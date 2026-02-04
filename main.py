import streamlit as st
import pandas as pd
import torch
import requests
import numpy as np
from Bio import SeqIO
from io import StringIO
from transformers import AutoTokenizer, EsmModel
import py3Dmol
from stmol import showmol

# --- ページ設定 ---
st.set_page_config(page_title="Vector2Fold - Protein Predictor", layout="wide", page_icon="🧬")

# --- モデルのキャッシュ ---
@st.cache_resource
def load_models():
    # CPU環境での動作を想定（GPUがある場合は "cuda" に変更可能）
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model_name = "facebook/esm2_t6_8M_UR50D"
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = EsmModel.from_pretrained(model_name).to(device)
    return tokenizer, model, device

# --- 関数群 ---
def translate_dna(dna_sequence):
    from Bio.Seq import Seq
    # 簡易的に最初のストップコドンまでを翻訳
    return str(Seq(dna_sequence).translate(to_stop=True))

def predict_solubility(protein_seq, tokenizer, model, device):
    # 特徴抽出による簡易スコアリング
    inputs = tokenizer(protein_seq, return_tensors="pt", padding=True, truncation=True).to(device)
    with torch.no_grad():
        outputs = model(**inputs)
        embeddings = outputs.last_hidden_state.mean(dim=1).cpu().numpy()
    # デモ用の擬似スコア（実際は学習済みモデルが必要）
    score = np.tanh(embeddings.mean()) * 50 + 50
    return float(score)

def get_esmfold_data(protein_seq):
    """
    ESMFold APIからPDBデータと平均pLDDTを取得する。
    pLDDTが0.0-1.0スケールの場合は自動的に0-100スケールへ正規化する。
    """
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    try:
        # 配列長制限（ESMFold APIの制限に合わせる）
        if len(protein_seq) > 1000:
             return None, None, "Sequence too long for API (>1000 AA)"

        response = requests.post(url, data=protein_seq, timeout=120)
        
        if response.status_code == 200:
            pdb_text = response.text
            # PDBのB-factor列(60-66文字目)から各残基のpLDDTを取得
            plddt_values = [float(line[60:66].strip()) for line in pdb_text.splitlines() if line.startswith("ATOM")]
            
            if not plddt_values:
                return pdb_text, None, "No pLDDT data found in PDB"

            # --- 自動正規化ロジックの開始 ---
            # 最大値が 1.0 以下の場合は、APIが小数スケールで返していると判断
            max_plddt = max(plddt_values)
            if max_plddt <= 1.0:
                # すべての値を100倍して 0-100 の範囲に変換
                plddt_values = [v * 100 for v in plddt_values]
            # --- 自動正規化ロジックの終了 ---

            # 平均値を計算
            avg_plddt = sum(plddt_values) / len(plddt_values)
            
            return pdb_text, avg_plddt, "Success"
        else:
            return None, None, f"API Error: {response.status_code}"
            
    except Exception as e:
        return None, None, str(e)
    
def render_mol(pdb_data):
    """py3Dmolを使って構造を表示する設定"""
    view = py3Dmol.view(width=800, height=600)
    view.addModel(pdb_data, 'pdb')
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    # pLDDTに基づいた色付け（青=高信頼度(90), 赤=低信頼度(50)）
    view.setStyle({'cartoon': {'colorscheme': {'prop':'b','gradient': 'roygb','min':0.5, 'max':0.9}}})
    view.zoomTo()
    return view

# --- UIセクション ---
st.title("🧬 Vector2Fold: 構造予測 & 可視化")
st.markdown("""
DNA配列からタンパク質の溶解性とフォールディングを予測します。
ESMFoldによって予測された3D構造は、pLDDTスコア（信頼度）に基づいて色分け表示されます。
(青: 高信頼度 > 緑 > 黄 > 赤: 低信頼度/ディスオーダー領域)
""")

# サイドバー
st.sidebar.header("Input Settings")
input_mode = st.sidebar.radio("入力方法:", ["テキスト入力", "FASTAファイルアップロード"])

dna_input = ""
if input_mode == "テキスト入力":
    dna_input = st.text_area("DNA塩基配列を入力 (複数可, FASTA形式):", height=200, placeholder=">Seq1\nATGGCC...")
else:
    uploaded_file = st.sidebar.file_uploader("FASTAファイルを選択", type=["fasta", "fa"])
    if uploaded_file:
        dna_input = uploaded_file.getvalue().decode("utf-8")

# 解析実行
if st.button("解析開始", type="primary"):
    if not dna_input or not ">" in dna_input:
        st.error("有効なFASTA形式の配列を入力してください。")
    else:
        # モデルロード（キャッシュ使用）
        with st.spinner("モデルを準備中..."):
            tokenizer, model, device = load_models()
        
        results = []
        fasta_io = StringIO(dna_input)
        records = list(SeqIO.parse(fasta_io, "fasta"))
        
        if not records:
             st.error("配列が見つかりませんでした。")
        else:
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            # 各配列の解析ループ
            for i, record in enumerate(records):
                status_text.text(f"解析中 ({i+1}/{len(records)}): {record.id}...")
                protein_seq = translate_dna(str(record.seq))
                
                # 予測実行
                sol_score = predict_solubility(protein_seq, tokenizer, model, device)
                pdb_data, plddt, api_status = get_esmfold_data(protein_seq)
                
                verdict = "✅ Confident" if plddt and plddt >= 70 else "⚠️ Caution/Disordered"
                if not plddt: verdict = f"❌ Failed ({api_status})"
                
                results.append({
                    "ID": record.id,
                    "Length": len(protein_seq),
                    "Solubility": f"{sol_score:.1f}",
                    "pLDDT": f"{plddt:.1f}" if plddt else "N/A",
                    "Verdict": verdict,
                    "PDB_Data": pdb_data # 構造データを保持
                })
                progress_bar.progress((i + 1) / len(records))
            
            status_text.empty()
            progress_bar.empty()
            st.success("解析完了！")

            # 結果のDataFrame化と表示
            df = pd.DataFrame(results)
            st.subheader("📊 解析結果一覧")
            st.dataframe(df.drop(columns=["PDB_Data"]), use_container_width=True) # PDBデータはテーブルには表示しない

            st.divider()

            # --- 3D構造ビューワーセクション ---
            st.subheader("🧊 3D構造ビューワー")
            
            # PDBデータ取得に成功した配列のみを選択肢にする
            available_structures = df[df["PDB_Data"].notnull()]
            
            if not available_structures.empty:
                selected_id = st.selectbox(
                    "表示するタンパク質を選択してください:",
                    available_structures["ID"].tolist()
                )
                
                # 選択されたIDに対応するPDBデータを取得
                selected_pdb = available_structures[available_structures["ID"] == selected_id]["PDB_Data"].iloc[0]
                selected_plddt = available_structures[available_structures["ID"] == selected_id]["pLDDT"].iloc[0]

                col1, col2 = st.columns([3, 1])
                with col1:
                    # py3Dmolによるレンダリング表示
                    view = render_mol(selected_pdb)
                    showmol(view, height=500)
                    st.caption("マウス操作: 回転(左ドラッグ), ズーム(ホイール), 移動(右ドラッグ)")
                with col2:
                    st.metric("選択中の配列 pLDDT", selected_plddt)
                    st.info("青色が濃いほど予測の信頼度が高く、安定した構造を形成する可能性が高い領域です。赤色はディスオーダー領域や柔軟性が高い領域を示唆します。")
                    # PDBダウンロードボタン
                    st.download_button(
                        label="PDBファイルをダウンロード",
                        data=selected_pdb,
                        file_name=f"{selected_id}_predicted.pdb",
                        mime="chemical/x-pdb"
                    )
            else:
                st.warning("表示可能な3D構造データがありません。")