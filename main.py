import streamlit as st
import pandas as pd
import requests
import io
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import Analysis, AllEnzymes
from dna_features_viewer import BiopythonTranslator
import py3Dmol
from stmol import showmol

# --- 既存の解析関数 ---

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
    view.setStyle({'cartoon': {'colorscheme': {'prop':'b', 'gradient': 'roygb', 'min': 50, 'max': 90}}})
    view.zoomTo()
    return view

# --- 新規：ベクター設計・出力用関数 ---

def generate_vector_map(sequence, label="New Vector"):
    record = SeqRecord(Seq(sequence), id="Vector", name=label)
    # 仮の特徴量（Insert部分）を目立たせるための設定
    translator = BiopythonTranslator()
    graphic_record = translator.translate_record(record)
    fig, ax = plt.subplots(figsize=(8, 8))
    graphic_record.plot(ax=ax, with_ruler=True)
    return fig

def export_to_excel(dna_seq, protein_seq, plddt, map_fig):
    output = io.BytesIO()
    map_img = io.BytesIO()
    map_fig.savefig(map_img, format='png', bbox_inches='tight')
    
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        df = pd.DataFrame({
            "項目": ["DNA配列長", "タンパク質配列", "予測平均pLDDT"],
            "値": [len(dna_seq), protein_seq, plddt]
        })
        df.to_excel(writer, sheet_name='解析レポート', index=False)
        
        # 配列詳細シート
        df_seq = pd.DataFrame({"DNA配列": [dna_seq]})
        df_seq.to_excel(writer, sheet_name='配列詳細', index=False)
        
        # 画像の挿入
        worksheet = writer.sheets['解析レポート']
        worksheet.insert_image('D2', 'map.png', {'image_data': map_img})
        
    return output.getvalue()

# --- Streamlit UI ---

st.set_page_config(page_title="Vector2Fold - Integrated Pipeline", layout="wide")
st.title("🧬 Vector2Fold: ベクター設計 & 構造予測パイプライン")

# サイドバー：ファイルアップロード
st.sidebar.header("1. 配列データの読み込み")
gene_file = st.sidebar.file_uploader("目的遺伝子 (Insert) FASTA/GB", type=["fasta", "gb", "txt"])
vector_file = st.sidebar.file_uploader("ベクター (Backbone) FASTA/GB", type=["fasta", "gb", "txt"])

if gene_file and vector_file:
    # 配列の読み込み
    gene_rec = SeqIO.read(io.StringIO(gene_file.getvalue().decode("utf-8")), "fasta")
    vector_rec = SeqIO.read(io.StringIO(vector_file.getvalue().decode("utf-8")), "fasta")
    
    st.subheader("Step 1: インシリコ・クローニング")
    
    col1, col2 = st.columns(2)
    with col1:
        st.info(f"目的遺伝子: {gene_rec.id} ({len(gene_rec.seq)} bp)")
        st.info(f"バックボーン: {vector_rec.id} ({len(vector_rec.seq)} bp)")
        
        # 制限酵素解析
        analysis = Analysis(AllEnzymes, vector_rec.seq)
        unique_sites = list(analysis.unique_sites().keys())
        selected_site = st.selectbox("挿入に使用するユニーク制限酵素サイトの候補:", unique_sites)
        
        insert_pos = st.number_input("挿入開始位置 (bp) ※手動調整", value=0, max_value=len(vector_rec.seq))
        
    # クローニング後の配列生成
    final_dna_seq = vector_rec.seq[:insert_pos] + gene_rec.seq + vector_rec.seq[insert_pos:]
    
    with col2:
        st.write("新ベクターマップ（簡易表示）")
        fig = generate_vector_map(str(final_dna_seq))
        st.pyplot(fig)

    st.divider()
    st.subheader("Step 2: タンパク質翻訳 & 構造予測")
    
    # 翻訳（目的遺伝子部分を想定して翻訳）
    protein_seq = str(gene_rec.seq.translate(to_stop=True))
    st.text_area("翻訳されたアミノ酸配列", protein_seq, height=100)
    
    if st.button("構造解析を開始"):
        with st.spinner("ESMFoldで構造を予測中..."):
            pdb_data, plddt, msg = get_esmfold_data(protein_seq)
            
            if pdb_data:
                res_col1, res_col2 = st.columns([1, 2])
                with res_col1:
                    st.metric("平均 pLDDT", f"{plddt:.2f}")
                    if plddt > 70:
                        st.success("高い信頼度で構造が予測されました。")
                    else:
                        st.warning("信頼度が低いため、天然変性領域の可能性があります。")
                    
                    # Excelダウンロード
                    excel_data = export_to_excel(str(final_dna_seq), protein_seq, plddt, fig)
                    st.download_button(
                        label="解析結果をExcelでダウンロード",
                        data=excel_data,
                        file_name="Vector2Fold_Report.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                    )
                
                with res_col2:
                    st.write("予測3D構造 (Color by pLDDT)")
                    showmol(render_mol(pdb_data), height=500, width=700)
            else:
                st.error(f"解析に失敗しました: {msg}")

else:
    st.write("サイドバーから目的遺伝子とベクターの配列ファイルをアップロードしてください。")