"""
Vector2Fold: DNA Sequence to Protein 3D Structure Pipeline
ベクターDNA配列 → タンパク質翻訳 → 物性解析 → 3D構造予測 → ベンチプロトコル生成
"""

import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import requests
import io
import math
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Restriction import Analysis, AllEnzymes, RestrictionBatch
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from dna_features_viewer import BiopythonTranslator, CircularGraphicRecord
import py3Dmol

from config import (
    MAJOR_ENZYMES,
    ESMFOLD_URL,
    ESMFOLD_MAX_LENGTH,
    ESMFOLD_TIMEOUT,
    PCR_PROTOCOLS,
    ENZYME_CONDITIONS,
)
from protocol_generator import generate_protocol_pdf, calc_tm_simple


# =========================================================================
# 解析ロジック
# =========================================================================

def get_esmfold_data(protein_seq: str) -> tuple:
    """ESMFold APIによる構造予測（3段階フォールバック）"""
    if len(protein_seq) > ESMFOLD_MAX_LENGTH:
        return None, None, f"Sequence too long (>{ESMFOLD_MAX_LENGTH} AA)"

    # Stage 1: verify=True
    for verify in [True, False]:
        try:
            response = requests.post(
                ESMFOLD_URL, data=protein_seq,
                timeout=ESMFOLD_TIMEOUT, verify=verify,
            )
            if response.status_code == 200:
                pdb_text = response.text
                plddt_values = [
                    float(line[60:66].strip())
                    for line in pdb_text.splitlines()
                    if line.startswith("ATOM")
                ]
                avg_plddt = (
                    (sum(plddt_values) / len(plddt_values)) * 100
                    if plddt_values else None
                )
                return pdb_text, avg_plddt, "Success"
            if verify:
                continue
            return None, None, f"API Error: {response.status_code}"
        except requests.exceptions.SSLError:
            if verify:
                continue
            return None, None, "SSL Error"
        except Exception as e:
            return None, None, str(e)

    # Stage 3: AlphaFold EBI fallback (UniProt IDが必要なため、ここではスキップ)
    return None, None, "ESMFold API unavailable"


def calculate_protein_properties(protein_seq: str) -> tuple:
    """アミノ酸配列から分子量・等電点・アミノ酸組成を計算"""
    cleaned = protein_seq.replace("*", "").replace("X", "")
    analysed = ProteinAnalysis(cleaned)
    mw = analysed.molecular_weight()
    pi = analysed.isoelectric_point()
    gravy = analysed.gravy()
    instability = analysed.instability_index()
    return mw, pi, gravy, instability


def render_mol_html(pdb_data: str, width: int = 700, height: int = 500) -> str:
    """py3Dmolの3D表示HTMLを生成（stmol不使用・直接HTML出力）"""
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_data, "pdb")
    view.setStyle({
        "cartoon": {
            "colorscheme": {
                "prop": "b",
                "gradient": "roygb",
                "min": 0.5,
                "max": 0.9,
            }
        }
    })
    view.zoomTo()
    return view._make_html()


def generate_vector_map(
    final_seq: str,
    insert_pos: int,
    insert_len: int,
    analysis_results=None,
    label: str = "New_Construct",
):
    """サーキュラーベクターマップを生成"""
    record = SeqRecord(Seq(final_seq), id="Vector", name=label)
    features = [
        SeqFeature(
            FeatureLocation(0, insert_pos),
            type="backbone",
            qualifiers={"label": "Backbone_5'", "color": "#e8e8e8"},
        ),
        SeqFeature(
            FeatureLocation(insert_pos, insert_pos + insert_len),
            type="insert",
            qualifiers={"label": "TARGET GENE", "color": "#ff8c00"},
        ),
        SeqFeature(
            FeatureLocation(insert_pos + insert_len, len(final_seq)),
            type="backbone",
            qualifiers={"label": "Backbone_3'", "color": "#e8e8e8"},
        ),
    ]

    if analysis_results:
        for enzyme, cuts in analysis_results.items():
            if enzyme in MAJOR_ENZYMES and len(cuts) == 1:
                pos = cuts[0]
                features.append(
                    SeqFeature(
                        FeatureLocation(pos - 1, pos),
                        type="restriction_site",
                        qualifiers={"label": str(enzyme), "color": "#d62728"},
                    )
                )

    record.features = features
    translator = BiopythonTranslator()
    graphic_record = translator.translate_record(record)
    circular_rec = CircularGraphicRecord(
        sequence_length=len(final_seq),
        features=graphic_record.features,
    )
    fig, ax = plt.subplots(figsize=(8, 8))
    circular_rec.plot(ax=ax)
    plt.tight_layout()
    return fig


def export_to_excel(
    dna_seq, protein_seq, plddt, mw, pi, gravy, instability, map_fig
):
    """解析結果をExcelファイルとしてエクスポート"""
    output = io.BytesIO()
    map_img = io.BytesIO()
    map_fig.savefig(map_img, format="png", bbox_inches="tight", dpi=150)

    with pd.ExcelWriter(output, engine="xlsxwriter") as writer:
        df_summary = pd.DataFrame({
            "Item": [
                "DNA length (bp)",
                "Protein length (AA)",
                "Avg pLDDT",
                "MW (kDa)",
                "pI",
                "GRAVY index",
                "Instability index",
            ],
            "Value": [
                len(dna_seq),
                len(protein_seq),
                f"{plddt:.2f}" if plddt else "N/A",
                f"{mw / 1000:.2f}",
                f"{pi:.2f}",
                f"{gravy:.3f}",
                f"{instability:.2f}",
            ],
        })
        df_summary.to_excel(writer, sheet_name="Summary", index=False)

        ws = writer.sheets["Summary"]
        ws.insert_image("D2", "map.png", {"image_data": map_img})

        df_seq = pd.DataFrame({
            "DNA_Sequence": [dna_seq],
            "Protein_Sequence": [protein_seq],
        })
        df_seq.to_excel(writer, sheet_name="Sequences", index=False)

    return output.getvalue()


# =========================================================================
# Streamlit UI
# =========================================================================

st.set_page_config(
    page_title="Vector2Fold",
    page_icon="🧬",
    layout="wide",
)
st.title("🧬 Vector2Fold: DNA → Protein Structure & Bench Protocol")

# --- サイドバー: ファイルアップロード ---
st.sidebar.header("1. Upload Sequences")
gene_file = st.sidebar.file_uploader(
    "Insert gene (FASTA/GenBank)", type=["fasta", "gb", "txt"]
)
vector_file = st.sidebar.file_uploader(
    "Vector backbone (FASTA/GenBank)", type=["fasta", "gb", "txt"]
)

st.sidebar.divider()
st.sidebar.header("2. Options")
pcr_protocol = st.sidebar.selectbox(
    "PCR Protocol",
    list(PCR_PROTOCOLS.keys()),
    format_func=lambda k: PCR_PROTOCOLS[k]["name"],
)

primer_f = st.sidebar.text_input("Forward primer (optional)", placeholder="e.g. ATGCGATCG...")
primer_r = st.sidebar.text_input("Reverse primer (optional)", placeholder="e.g. TTACGATCG...")

# --- メインコンテンツ ---
if gene_file and vector_file:
    try:
        # 配列読み込み
        gene_content = gene_file.getvalue().decode("utf-8")
        vector_content = vector_file.getvalue().decode("utf-8")

        fmt = "genbank" if gene_file.name.endswith(".gb") else "fasta"
        gene_rec = SeqIO.read(io.StringIO(gene_content), fmt)

        fmt_v = "genbank" if vector_file.name.endswith(".gb") else "fasta"
        vector_rec = SeqIO.read(io.StringIO(vector_content), fmt_v)

        # ===== タブ構成 =====
        tab_clone, tab_analysis, tab_protocol = st.tabs([
            "🔬 In-silico Cloning",
            "🧪 Protein Analysis & 3D",
            "📋 Bench Protocol",
        ])

        # ===== Tab 1: クローニング =====
        with tab_clone:
            st.subheader("Step 1: In-silico Cloning")
            col1, col2 = st.columns([1, 2])

            with col1:
                analysis = Analysis(AllEnzymes, vector_rec.seq)
                all_cuts = analysis.full()
                unique_sites = sorted([
                    str(enz) for enz, cuts in all_cuts.items()
                    if len(cuts) == 1
                ])

                selected_site = None
                insert_pos = 0

                if unique_sites:
                    selected_site = st.selectbox(
                        "Unique restriction sites on vector:", unique_sites
                    )
                    target_enzyme = [
                        e for e in AllEnzymes if str(e) == selected_site
                    ][0]
                    suggested_pos = target_enzyme.search(vector_rec.seq)[0]
                    insert_pos = st.number_input(
                        "Insertion position (bp)",
                        value=suggested_pos,
                        max_value=len(vector_rec.seq),
                    )

                    # 制限酵素情報の表示
                    enz_info = ENZYME_CONDITIONS.get(selected_site, {})
                    if enz_info:
                        st.caption(f"Overhang: {enz_info['overhang']}")
                        st.caption(f"Buffer: {enz_info['buffer']}")
                        st.caption(f"Star activity: {enz_info['star_risk']}")
                else:
                    insert_pos = st.number_input("Insertion position (bp)", value=0)

                st.metric("Insert length", f"{len(gene_rec.seq)} bp")
                st.metric("Vector length", f"{len(vector_rec.seq)} bp")

            # コンストラクト構築
            final_dna_seq = (
                str(vector_rec.seq[:insert_pos])
                + str(gene_rec.seq)
                + str(vector_rec.seq[insert_pos:])
            )

            with col2:
                fig_map = generate_vector_map(
                    final_dna_seq, insert_pos, len(gene_rec.seq), all_cuts
                )
                st.pyplot(fig_map)

        # ===== Tab 2: タンパク質解析 =====
        with tab_analysis:
            st.subheader("Step 2: Protein Analysis & Structure Prediction")

            protein_seq = str(gene_rec.seq.translate(to_stop=True))
            st.text_area("Translated amino acid sequence", protein_seq, height=100)

            if st.button("Run Analysis", type="primary"):
                with st.spinner("Computing properties & predicting structure..."):
                    # 物性計算
                    mw, pi, gravy, instability = calculate_protein_properties(protein_seq)
                    # 構造予測
                    pdb_data, plddt, msg = get_esmfold_data(protein_seq)

                    # セッションに保存
                    st.session_state["analysis_done"] = True
                    st.session_state["mw"] = mw
                    st.session_state["pi"] = pi
                    st.session_state["gravy"] = gravy
                    st.session_state["instability"] = instability
                    st.session_state["pdb_data"] = pdb_data
                    st.session_state["plddt"] = plddt
                    st.session_state["protein_seq"] = protein_seq
                    st.session_state["final_dna_seq"] = final_dna_seq
                    st.session_state["fig_map"] = fig_map
                    st.session_state["msg"] = msg

            # 結果表示
            if st.session_state.get("analysis_done"):
                mw = st.session_state["mw"]
                pi = st.session_state["pi"]
                gravy = st.session_state["gravy"]
                instability = st.session_state["instability"]
                pdb_data = st.session_state["pdb_data"]
                plddt = st.session_state["plddt"]

                m1, m2, m3, m4 = st.columns(4)
                m1.metric("Avg pLDDT", f"{plddt:.1f}" if plddt else "N/A")
                m2.metric("MW", f"{mw / 1000:.2f} kDa")
                m3.metric("pI", f"{pi:.2f}")
                m4.metric("GRAVY", f"{gravy:.3f}")

                stability_label = "Stable" if instability < 40 else "Potentially unstable"
                st.caption(f"Instability index: {instability:.2f} ({stability_label})")

                if pdb_data:
                    res_col1, res_col2 = st.columns([1, 2])

                    with res_col1:
                        excel_data = export_to_excel(
                            final_dna_seq, protein_seq,
                            plddt, mw, pi, gravy, instability, fig_map,
                        )
                        st.download_button(
                            label="Download Excel Report",
                            data=excel_data,
                            file_name="Vector2Fold_Report.xlsx",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        )

                    with res_col2:
                        st.write("Predicted 3D Structure (colored by pLDDT)")
                        mol_html = render_mol_html(pdb_data)
                        components.html(mol_html, height=520, width=720, scrolling=False)
                else:
                    st.error(f"Structure prediction error: {st.session_state['msg']}")

        # ===== Tab 3: ベンチプロトコル =====
        with tab_protocol:
            st.subheader("Step 3: Generate Bench Protocol")

            st.write(
                "クローニング結果と解析データに基づいて、"
                "実験ベンチ用のプロトコルPDFを自動生成します。"
            )

            proto_col1, proto_col2 = st.columns(2)
            with proto_col1:
                construct_name = st.text_input(
                    "Construct name", value=f"{gene_rec.id}_in_{vector_rec.id}"
                )
            with proto_col2:
                st.write("")  # spacer
                st.caption(f"PCR protocol: {PCR_PROTOCOLS[pcr_protocol]['name']}")
                if primer_f:
                    tm_f = calc_tm_simple(primer_f)
                    tm_r = calc_tm_simple(primer_r) if primer_r else tm_f
                    st.caption(f"Estimated Tm: Fwd={tm_f:.0f} C / Rev={tm_r:.0f} C")

            if st.button("Generate Protocol PDF", type="primary"):
                with st.spinner("Generating protocol..."):
                    # セッション情報の取得
                    _mw = st.session_state.get("mw")
                    _pi = st.session_state.get("pi")
                    _plddt = st.session_state.get("plddt")
                    _protein = st.session_state.get("protein_seq", str(gene_rec.seq.translate(to_stop=True)))
                    _fig = st.session_state.get("fig_map", fig_map)

                    pdf_bytes = generate_protocol_pdf(
                        construct_name=construct_name,
                        insert_name=gene_rec.id,
                        vector_name=vector_rec.id,
                        insert_len_bp=len(gene_rec.seq),
                        vector_len_bp=len(vector_rec.seq),
                        selected_enzyme=selected_site,
                        insert_pos=insert_pos,
                        protein_seq=_protein,
                        mw_kda=_mw / 1000 if _mw else None,
                        pi_value=_pi,
                        plddt=_plddt,
                        primer_f=primer_f if primer_f else None,
                        primer_r=primer_r if primer_r else None,
                        protocol_key=pcr_protocol,
                        vector_map_fig=_fig,
                    )

                    st.download_button(
                        label="Download Protocol PDF",
                        data=pdf_bytes,
                        file_name=f"{construct_name}_protocol.pdf",
                        mime="application/pdf",
                        type="primary",
                    )

                    st.success("Protocol PDF generated successfully!")

                    # プロトコルのプレビュー表示
                    st.divider()
                    st.write("**Protocol Preview**")

                    preview_cols = st.columns(3)
                    with preview_cols[0]:
                        st.markdown("**PCR Conditions**")
                        proto_info = PCR_PROTOCOLS[pcr_protocol]
                        insert_kb = len(gene_rec.seq) / 1000
                        ext_fast = math.ceil(insert_kb * proto_info["extension_rate_fast"])
                        st.code(
                            f"{proto_info['denature_temp']}C  {proto_info['denature_time']} sec\n"
                            f"{proto_info['anneal_temp']}C  5-15 sec\n"
                            f"{proto_info['extension_temp']}C  {max(ext_fast,5)} sec\n"
                            f"x {proto_info['cycles_min']}-{proto_info['cycles_max']} cycles",
                            language=None,
                        )

                    with preview_cols[1]:
                        if selected_site:
                            st.markdown("**Restriction Enzyme**")
                            enz = ENZYME_CONDITIONS.get(selected_site, {})
                            st.code(
                                f"Enzyme: {selected_site}\n"
                                f"Overhang: {enz.get('overhang', '-')}\n"
                                f"Buffer: {enz.get('buffer', '-')}\n"
                                f"Temp: {enz.get('temp', 37)}C",
                                language=None,
                            )

                    with preview_cols[2]:
                        if _protein:
                            st.markdown("**Protein Summary**")
                            st.code(
                                f"Length: {len(_protein)} AA\n"
                                f"MW: {_mw/1000:.1f} kDa\n" if _mw else ""
                                f"pI: {_pi:.2f}\n" if _pi else ""
                                f"pLDDT: {_plddt:.1f}" if _plddt else "pLDDT: N/A",
                                language=None,
                            )

    except Exception as e:
        st.error(f"Error: {e}")
        st.exception(e)
else:
    st.info("👈 Upload insert gene and vector FASTA files from the sidebar to begin.")
    st.divider()

    st.subheader("About Vector2Fold")
    st.write(
        "Vector2Fold is a web application for researchers that integrates "
        "in-silico cloning, protein property analysis, ESMFold-based 3D "
        "structure prediction, and automatic bench protocol generation "
        "into a single pipeline."
    )

    feat_cols = st.columns(3)
    with feat_cols[0]:
        st.markdown("**🔬 In-silico Cloning**")
        st.caption("Vector map visualization with restriction enzyme analysis")
    with feat_cols[1]:
        st.markdown("**🧪 Protein Analysis**")
        st.caption("MW, pI, GRAVY, instability index, ESMFold 3D prediction")
    with feat_cols[2]:
        st.markdown("**📋 Bench Protocol**")
        st.caption("Auto-generated PDF with PCR, digestion, ligation, transformation")
