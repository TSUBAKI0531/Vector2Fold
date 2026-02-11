import streamlit as st
import pandas as pd
import requests
import io
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Restriction import Analysis, AllEnzymes, RestrictionBatch, EcoRI, BamHI, HindIII, NotI, XhoI, SpeI, PstI, NcoI, NdeI, SalI, KpnI, SacI, SmaI
# 追加: タンパク質の物性解析用モジュール
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from dna_features_viewer import BiopythonTranslator, CircularGraphicRecord
import py3Dmol
from stmol import showmol

# 主要な制限酵素のリスト
MajorEnzymes = RestrictionBatch([EcoRI, BamHI, HindIII, NotI, XhoI, SpeI, PstI, NcoI, NdeI, SalI, KpnI, SacI, SmaI])

# --- 解析ロジック ---

def get_esmfold_data(protein_seq):
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    try:
        if len(protein_seq) > 1000:
            return None, None, "Sequence too long (>1000 AA)"
        response = requests.post(url, data=protein_seq, timeout=120)
        if response.status_code == 200:
            pdb_text = response.text
            plddt_values = [float(line[60:66].strip()) for line in pdb_text.splitlines() if line.startswith("ATOM")]
            avg_plddt = (sum(plddt_values) / len(plddt_values)) * 100 if plddt_values else None
            return pdb_text, avg_plddt, "Success"
        return None, None, f"API Error: {response.status_code}"
    except Exception as e:
        return None, None, str(e)

def calculate_protein_properties(protein_seq):
    """アミノ酸配列から分子量と等電点を計算する"""
    analysed_seq = ProteinAnalysis(protein_seq)
    mw = analysed_seq.molecular_weight()
    pi = analysed_seq.isoelectric_point()
    return mw, pi

def render_mol(pdb_data):
    view = py3Dmol.view(width=800, height=600)
    view.addModel(pdb_data, 'pdb')
    # スケール修正: 0.5(赤)〜0.9(青)
    view.setStyle({'cartoon': {'colorscheme': {'prop':'b', 'gradient': 'roygb', 'min': 0.5, 'max': 0.9}}})
    view.zoomTo()
    return view

def generate_vector_map(final_seq, insert_pos, insert_len, analysis_results=None, label="New_Construct"):
    record = SeqRecord(Seq(final_seq), id="Vector", name=label)
    features = []
    features.append(SeqFeature(FeatureLocation(0, insert_pos), type="backbone", qualifiers={"label": "Backbone_A", "color": "#f4f4f4"}))
    features.append(SeqFeature(FeatureLocation(insert_pos, insert_pos + insert_len), type="insert", qualifiers={"label": "TARGET GENE", "color": "#ffaa00"}))
    features.append(SeqFeature(FeatureLocation(insert_pos + insert_len, len(final_seq)), type="backbone", qualifiers={"label": "Backbone_B", "color": "#f4f4f4"}))
    
    if analysis_results:
        for enzyme, cuts in analysis_results.items():
            if enzyme in MajorEnzymes and len(cuts) == 1:
                pos = cuts[0]
                actual_pos = pos if pos <= insert_pos else pos + insert_len
                features.append(SeqFeature(FeatureLocation(actual_pos-1, actual_pos), type="restriction_site", qualifiers={"label": str(enzyme), "color": "#d62728"}))
    
    record.features = features
    translator = BiopythonTranslator()
    graphic_record = translator.translate_record(record)
    circular_rec = CircularGraphicRecord(sequence_length=len(final_seq), features=graphic_record.features)
    fig, ax = plt.subplots(figsize=(10, 10))
    circular_rec.plot(ax=ax)
    return fig

def export_to_excel(dna_seq, protein_seq, plddt, mw, pi, map_fig):
    output = io.BytesIO()
    map_img = io.BytesIO()
    map_fig.savefig(map_img, format='png', bbox_inches='tight')
    
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        df = pd.DataFrame({
            "項目": ["DNA配列長(bp)", "アミノ酸配列", "平均pLDDT", "分子量(kDa)", "等電点(pI)"],
            "値": [len(dna_seq), protein_seq, f"{plddt:.2f}", f"{mw/1000:.2f}", f"{pi:.2f}"]
        })
        df.to_excel(writer, sheet_name='Summary', index=False)
        worksheet = writer.sheets['Summary']
        worksheet.insert_image('D2', 'map.png', {'image_data': map_img})
        df_seq = pd.DataFrame({"DNA_Sequence": [dna_seq]})
        df_seq.to_excel(writer, sheet_name='Sequence_Detail', index=False)
        
    return output.getvalue()

# --- UI ---

st.set_page_config(page_title="Vector2Fold", layout="wide")
st.title("🧬 Vector2Fold: 研究用タンパク質解析ツール")

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
            analysis = Analysis(AllEnzymes, vector_rec.seq)
            all_cuts = analysis.full()
            unique_sites = sorted([str(enzyme) for enzyme, cuts in all_cuts.items() if len(cuts) == 1])
            
            if unique_sites:
                selected_site = st.selectbox("利用可能なユニーク制限酵素サイト:", unique_sites)
                target_enzyme = [e for e in AllEnzymes if str(e) == selected_site][0]
                suggested_pos = target_enzyme.search(vector_rec.seq)[0]
                insert_pos = st.number_input("挿入位置 (bp) の指定", value=suggested_pos, max_value=len(vector_rec.seq))
            else:
                insert_pos = st.number_input("挿入位置 (bp) の指定", value=0)
        
        final_dna_seq = str(vector_rec.seq[:insert_pos]) + str(gene_rec.seq) + str(vector_rec.seq[insert_pos:])
        
        with col2:
            fig = generate_vector_map(final_dna_seq, insert_pos, len(gene_rec.seq), all_cuts)
            st.pyplot(fig)

        st.divider()
        st.subheader("Step 2: タンパク質物性解析 & 構造予測")
        
        protein_seq = str(gene_rec.seq.translate(to_stop=True))
        st.text_area("翻訳後のアミノ酸配列", protein_seq, height=100)
        
        if st.button("総合解析を開始"):
            with st.spinner("物性計算および構造予測中..."):
                # 物性計算
                mw, pi = calculate_protein_properties(protein_seq)
                # 構造予測
                pdb_data, plddt, msg = get_esmfold_data(protein_seq)
                
                if pdb_data:
                    m1, m2, m3 = st.columns(3)
                    m1.metric("平均 pLDDT", f"{plddt:.2f}")
                    m2.metric("分子量 (MW)", f"{mw/1000:.2f} kDa")
                    m3.metric("等電点 (pI)", f"{pi:.2f}")
                    
                    res_col1, res_col2 = st.columns([1, 2])
                    with res_col1:
                        excel_data = export_to_excel(final_dna_seq, protein_seq, plddt, mw, pi, fig)
                        st.download_button(label="レポートをExcelでダウンロード", data=excel_data, file_name="Vector2Fold_Report.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
                    with res_col2:
                        st.write("予測3D構造 (Color by pLDDT)")
                        showmol(render_mol(pdb_data), height=500, width=700)
                else:
                    st.error(f"解析エラー: {msg}")
    except Exception as e:
        st.error(f"エラーが発生しました: {e}")
else:
    st.write("サイドバーからFASTAファイルをアップロードしてください。")