# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running the App

```bash
pip install -r requirements.txt
streamlit run main.py
```

## Architecture

Three-module design:

- **`main.py`** — Streamlit UI and orchestration. All tab layout, session state management, and calls to analysis/generation functions live here.
- **`config.py`** — All constants: restriction enzyme database (`ENZYME_CONDITIONS`), PCR protocol parameters (`PCR_PROTOCOLS`), ESMFold API settings, and `MAJOR_ENZYMES` batch. **Values must remain ASCII-only** (fpdf2 Helvetica Latin-1 constraint).
- **`protocol_generator.py`** — fpdf2-based multi-page PDF generation. Contains `ProtocolPDF` class, `_s()` Unicode sanitizer, and `generate_protocol_pdf()` entry point.

## Pipeline Flow

User uploads insert gene FASTA/GenBank + vector backbone → three tabs:

1. **In-silico Cloning**: Biopython `Analysis(AllEnzymes)` finds unique restriction sites on the vector; user selects enzyme and insertion position; `generate_vector_map()` renders a circular map via `dna_features_viewer`.
2. **Protein Analysis**: Gene sequence translated via Biopython; `ProteinAnalysis` computes MW/pI/GRAVY/instability; `get_esmfold_data()` calls ESMFold API (2-stage SSL fallback); `render_mol_html()` generates py3Dmol HTML embedded via `streamlit.components.v1.html()`. Results stored in `st.session_state` to persist across tabs.
3. **Bench Protocol**: `generate_protocol_pdf()` builds a multi-page PDF (PCR, digestion, ligation, transformation, colony screening, troubleshooting).

## Key Implementation Details

- **No stmol dependency**: 3D rendering uses `py3Dmol._make_html()` + `components.html()` directly for Streamlit Cloud compatibility.
- **ESMFold fallback**: `get_esmfold_data()` tries `verify=True` then `verify=False`; returns `(None, None, error_msg)` if both fail. Max sequence length is 1000 AA (`ESMFOLD_MAX_LENGTH`).
- **fpdf2 Unicode**: `_s()` in `protocol_generator.py` maps common Unicode characters to ASCII equivalents before any PDF cell/text call. Strings in `config.py` are already sanitized at definition time (see NOTE comments).
- **Session state pattern**: Analysis results (`mw`, `pi`, `plddt`, `pdb_data`, `fig_map`, etc.) are written to `st.session_state` after "Run Analysis" and read back in the Protocol tab to avoid re-computation.
