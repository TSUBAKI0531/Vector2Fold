import streamlit as st
import pandas as pd
import requests
import io
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Restriction import Analysis, AllEnzymes
# 修正：CircularGraphicRecord を追加インポート
from dna_features_viewer import BiopythonTranslator, CircularGraphicRecord
import py3Dmol
from stmol import showmol

# --- 予測・描画ロジック ---

def get_esmfold_data(protein_seq):
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    try:
        if len(protein_seq) > 1000:
            return None, None, "Sequence too long (>1000 AA)"
        response = requests.post(url, data=protein_seq, timeout=120)
        if response.status_code == 200:
            pdb_text = response.text
            plddt_values = [float(line[60:66].strip()) for line in pdb_text.splitlines() if line.startswith("ATOM")]
            if plddt_values and max(plddt_values) <= 1.0:
                plddt_values = [v * 100 for v in plddt_values]
            avg_plddt = sum(plddt_values) / len(plddt_values) if plddt_values else None
            return pdb_text, avg_plddt, "Success"
        return None, None, f"API Error: {response.status_code}"
    except Exception as e:
        return None, None, str(e)

def render_mol(pdb_data):
    view = py3Dmol.view(width=800, height=600)
    view.addModel(pdb_data, 'pdb')
    # pLDDTに基づいた色分け: 50(赤)〜90(青)
    view.setStyle({'cartoon': {'colorscheme': {'prop':'b', 'gradient': 'roygb', 'min': 50, 'max': 90}}})
    view.zoomTo()
    return view

def generate_vector_map(final_seq, insert_pos, insert_len, label="New_Construct"):
    """円形マップを正しく生成する"""
    record = SeqRecord(Seq(final_seq), id="Vector", name=label)
    
    # 特徴量（アノテーション）の定義
    feat1 = SeqFeature(FeatureLocation(0, insert_pos), type="backbone", qualifiers={"label": "Backbone_A", "color": "#f4f4f4"})
    feat2 = SeqFeature(FeatureLocation(insert_pos, insert_pos + insert_len), type="insert", qualifiers={"label": "TARGET GENE", "color": "#ffaa00"})
    feat3 = SeqFeature(FeatureLocation(insert_pos + insert_len, len(final_seq)), type="backbone", qualifiers={"label": "Backbone_B", "color": "#f4f4f4"})
    record.features = [feat1, feat2, feat3]
    
    # Translatorを使って標準的なGraphicRecordを一度作成し、特徴量を抽出
    translator = BiopythonTranslator()
    graphic_record = translator.translate_record(record)
    
    # 【重要】円形マップ専用のクラスに変換
    circular_rec = CircularGraphicRecord(
        sequence_length=len(final_seq),
        features=graphic_record.features
    )
    
    fig, ax = plt.subplots(figsize=(8, 8))
    # plot() を呼び出すと、CircularGraphicRecord なので円形に描画される
    circular_rec.plot(ax=ax)
    return fig

def export_to_excel(dna_seq, protein_seq, plddt, map_fig):
    output = io.BytesIO()
    map_img = io.BytesIO()
    map_fig.savefig(map_img, format='png', bbox_inches='tight')
    
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        df = pd.DataFrame({
            "項目": ["DNA配列長", "アミノ酸配列", "予測平均pLDDT"],
            "値": [len(dna_seq), protein_seq, plddt]
        })
        df.to_excel(writer, sheet_name='Summary', index=False)
        worksheet = writer.sheets['Summary']
        worksheet.insert_image('D2', 'map.png', {'image_data': map_img})
        df_seq = pd.DataFrame({"DNA_Sequence": [dna_seq]})
        df_seq.to_excel(writer, sheet_name='Sequence_Detail', index=False)
        
    return output.getvalue()

# --- アプリケーションUI ---

st.set_page_config(page_title="Vector2Fold", layout="wide")
st.title("🧬 Vector2Fold: ベクター設計 & 構造予測")

st.sidebar.header("1. 配列データの読み込み")
gene_file = st.sidebar.file_uploader("目的遺伝子 (Insert) FASTA/GB", type=["fasta", "gb", "txt"])
vector_file = st.sidebar.file_uploader("ベクター (Backbone) FASTA/GB", type=["fasta", "gb", "txt"])

if gene_file and vector_file:
    try:
        gene_rec = SeqIO.read(io.StringIO(gene_file.getvalue().decode("utf-8")), "fasta")
        vector_rec = SeqIO.read(io.StringIO(vector_file.getvalue().decode("utf-8")), "fasta")
        
        st.subheader("Step 1: インシリコ・クローニング")
        
        col1, col2 = st.columns(2)
        with col1:
            st.info(f"目的遺伝子: {gene_rec.id} ({len(gene_rec.seq)} bp)")
            st.info(f"バックボーン: {vector_rec.id} ({len(vector_rec.seq)} bp)")
            
            # 制限酵素解析
            analysis = Analysis(AllEnzymes, vector_rec.seq)
            all_cuts = analysis.full()
            unique_sites = [str(enzyme) for enzyme, cuts in all_cuts.items() if len(cuts) == 1]
            unique_sites.sort()
            
            if unique_sites:
                selected_site = st.selectbox("利用可能なユニーク制限酵素サイト:", unique_sites)
            else:
                st.warning("ユニーク制限酵素サイトが見つかりませんでした。")
            
            insert_pos = st.number_input("挿入位置 (bp) の指定", value=0, max_value=len(vector_rec.seq))
        
        # 配列の結合
        final_dna_seq = str(vector_rec.seq[:insert_pos]) + str(gene_rec.seq) + str(vector_rec.seq[insert_pos:])
        
        with col2:
            st.write("プレビュー: 合成後のベクターマップ")
            # 修正した関数を呼び出し
            fig = generate_vector_map(final_dna_seq, insert_pos, len(gene_rec.seq))
            st.pyplot(fig)

        st.divider()
        st.subheader("Step 2: 構造予測とレポート出力")
        
        protein_seq = str(gene_rec.seq.translate(to_stop=True))
        st.text_area("翻訳後のアミノ酸配列", protein_seq, height=100)
        
        if st.button("立体構造の解析を開始"):
            with st.spinner("解析中..."):
                pdb_data, plddt, msg = get_esmfold_data(protein_seq)
                
                if pdb_data:
                    res_col1, res_col2 = st.columns([1, 2])
                    with res_col1:
                        st.metric("平均 pLDDT", f"{plddt:.2f}")
                        excel_data = export_to_excel(final_dna_seq, protein_seq, plddt, fig)
                        st.download_button(
                            label="レポートをExcelでダウンロード",
                            data=excel_data,
                            file_name="Vector2Fold_Report.xlsx",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                        )
                    with res_col2:
                        st.write("予測3D構造 (Color by pLDDT)")
                        showmol(render_mol(pdb_data), height=500, width=700)
                else:
                    st.error(f"解析エラー: {msg}")
    except Exception as e:
        st.error(f"エラーが発生しました: {e}")
else:
    st.write("サイドバーからFASTAファイルをアップロードしてください。")