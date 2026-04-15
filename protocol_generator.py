"""
Vector2Fold Protocol Auto-Generation Module
Generates bench protocol PDFs for cloning experiments.
"""

import io
import math
import re
from datetime import datetime
from fpdf import FPDF

from config import (
    PCR_PROTOCOLS,
    ENZYME_CONDITIONS,
    LIGATION_PROTOCOLS,
)


# =========================================================================
# Text sanitizer for fpdf2 Helvetica font (Latin-1 only)
# =========================================================================

_UNICODE_REPLACEMENTS = {
    "\u301C": "-",   # 〜 (wave dash)
    "\uFF5E": "-",   # ～ (fullwidth tilde)
    "\u00D7": "x",   # × (multiplication sign)
    "\u00B5": "u",   # µ (micro sign)
    "\u03BC": "u",   # μ (Greek mu)
    "\u2212": "-",   # − (minus sign)
    "\u2013": "-",   # – (en dash)
    "\u2014": "-",   # — (em dash)
    "\u2018": "'",   # ' (left single quote)
    "\u2019": "'",   # ' (right single quote)
    "\u201C": '"',   # " (left double quote)
    "\u201D": '"',   # " (right double quote)
    "\u2026": "...", # … (ellipsis)
    "\u00B0": " ",   # ° (degree sign) → use "C" instead of "°C"
    "\u2192": "->",  # → (right arrow)
    "\u2190": "<-",  # ← (left arrow)
    "\u2265": ">=",  # ≥
    "\u2264": "<=",  # ≤
}


def sanitize_for_pdf(text: str) -> str:
    """
    Remove or replace characters outside fpdf2 Helvetica (Latin-1) range.
    Applies known substitutions first, then strips any remaining non-Latin-1
    characters to prevent FPDFUnicodeEncodingException.
    """
    if not isinstance(text, str):
        text = str(text)
    # Known replacements
    for char, replacement in _UNICODE_REPLACEMENTS.items():
        text = text.replace(char, replacement)
    # Strip any remaining characters outside Latin-1 (U+0000-U+00FF)
    text = "".join(c if ord(c) <= 0x00FF else "?" for c in text)
    return text


# =========================================================================
# Custom PDF class
# =========================================================================

class ProtocolPDF(FPDF):
    """Protocol PDF with built-in text sanitization."""

    def __init__(self):
        super().__init__()
        self.set_auto_page_break(auto=True, margin=20)

    # --- Override core text methods to auto-sanitize ---

    def cell(self, w=0, h=0, text="", border=0, new_x=None, new_y=None,
             align="", fill=False, center=False, link="", markdown=False):
        return super().cell(
            w, h, sanitize_for_pdf(text), border, new_x, new_y,
            align, fill, center, link, markdown,
        )

    def multi_cell(self, w=0, h=0, text="", border=0, align="",
                   fill=False, new_x=None, new_y=None, **kwargs):
        return super().multi_cell(
            w, h, sanitize_for_pdf(text), border, align,
            fill, new_x, new_y, **kwargs,
        )

    # --- Header / Footer (use super() to avoid double-sanitize) ---

    def header(self):
        self.set_font("Helvetica", "B", 10)
        self.set_text_color(100, 100, 100)
        self.cell(0, 6, "Vector2Fold - Auto-generated Bench Protocol",
                  align="R", new_x="LMARGIN", new_y="NEXT")
        self.set_draw_color(0, 102, 204)
        self.set_line_width(0.5)
        self.line(10, self.get_y(), 200, self.get_y())
        self.ln(4)

    def footer(self):
        self.set_y(-15)
        self.set_font("Helvetica", "I", 8)
        self.set_text_color(128, 128, 128)
        self.cell(0, 10, f"Page {self.page_no()}/{{nb}}", align="C")

    # --- Convenience methods ---

    def section_title(self, title: str, level: int = 1):
        if level == 1:
            self.set_font("Helvetica", "B", 14)
            self.set_text_color(0, 51, 102)
            self.set_fill_color(230, 240, 250)
            self.cell(0, 10, title, fill=True, new_x="LMARGIN", new_y="NEXT")
        elif level == 2:
            self.set_font("Helvetica", "B", 11)
            self.set_text_color(0, 80, 140)
            self.cell(0, 8, title, new_x="LMARGIN", new_y="NEXT")
        self.ln(2)

    def body_text(self, text: str):
        self.set_font("Helvetica", "", 10)
        self.set_text_color(30, 30, 30)
        self.multi_cell(0, 5, text)
        self.ln(2)

    def note_box(self, text: str, box_type: str = "note"):
        colors = {
            "note":    (255, 255, 220, 200, 180, 50),
            "warning": (255, 230, 230, 200, 80, 80),
            "tip":     (220, 245, 220, 60, 160, 60),
        }
        bg_r, bg_g, bg_b, txt_r, txt_g, txt_b = colors.get(
            box_type, colors["note"])
        labels = {"note": "[NOTE]", "warning": "[WARNING]", "tip": "[TIP]"}

        self.set_fill_color(bg_r, bg_g, bg_b)
        self.set_draw_color(txt_r, txt_g, txt_b)

        self.set_font("Helvetica", "B", 9)
        self.set_text_color(txt_r, txt_g, txt_b)
        self.cell(0, 6, f"  {labels.get(box_type, '[NOTE]')}",
                  fill=True, new_x="LMARGIN", new_y="NEXT")
        self.set_font("Helvetica", "", 9)
        self.set_text_color(60, 60, 60)
        self.multi_cell(0, 5, f"  {text}", fill=True)
        self.ln(3)

    def add_table(self, headers: list, rows: list, col_widths: list = None):
        if col_widths is None:
            col_widths = [190 / len(headers)] * len(headers)

        # Header row
        self.set_font("Helvetica", "B", 9)
        self.set_fill_color(0, 80, 140)
        self.set_text_color(255, 255, 255)
        for i, header in enumerate(headers):
            self.cell(col_widths[i], 7, header,
                      border=1, fill=True, align="C")
        self.ln()

        # Data rows
        self.set_font("Helvetica", "", 9)
        self.set_text_color(30, 30, 30)
        for row_idx, row in enumerate(rows):
            if row_idx % 2 == 0:
                self.set_fill_color(245, 248, 252)
            else:
                self.set_fill_color(255, 255, 255)
            for i, cell_val in enumerate(row):
                self.cell(col_widths[i], 6, str(cell_val),
                          border=1, fill=True, align="C")
            self.ln()
        self.ln(3)


# =========================================================================
# Helper functions
# =========================================================================

def calc_tm_simple(primer_seq: str) -> float:
    """Simple Tm calculation: Tm = 2(A+T) + 4(G+C) - 5"""
    p = primer_seq.upper().strip()
    return (2 * (p.count("A") + p.count("T"))
            + 4 * (p.count("G") + p.count("C")) - 5)


def calc_pcr_extension_time(
    insert_len_bp: int, protocol_key: str = "PrimeSTAR_Max"
) -> dict:
    proto = PCR_PROTOCOLS[protocol_key]
    insert_kb = insert_len_bp / 1000
    fast_sec = math.ceil(insert_kb * proto["extension_rate_fast"])
    slow_sec = math.ceil(insert_kb * proto["extension_rate_slow"])
    return {
        "fast_sec": max(fast_sec, 5),
        "slow_sec": max(slow_sec, 10),
        "insert_kb": round(insert_kb, 2),
    }


# =========================================================================
# Main PDF generation
# =========================================================================

def generate_protocol_pdf(
    construct_name: str,
    insert_name: str,
    vector_name: str,
    insert_len_bp: int,
    vector_len_bp: int,
    selected_enzyme: str = None,
    insert_pos: int = None,
    protein_seq: str = None,
    mw_kda: float = None,
    pi_value: float = None,
    plddt: float = None,
    primer_f: str = None,
    primer_r: str = None,
    protocol_key: str = "PrimeSTAR_Max",
    vector_map_fig=None,
) -> bytes:
    """Generate a bench protocol PDF."""

    pdf = ProtocolPDF()
    pdf.alias_nb_pages()
    pdf.add_page()

    proto = PCR_PROTOCOLS[protocol_key]
    ext_times = calc_pcr_extension_time(insert_len_bp, protocol_key)
    total_len = insert_len_bp + vector_len_bp

    # ===== Title =====
    pdf.set_font("Helvetica", "B", 18)
    pdf.set_text_color(0, 51, 102)
    pdf.cell(0, 12, "Bench Protocol",
             align="C", new_x="LMARGIN", new_y="NEXT")
    pdf.set_font("Helvetica", "", 11)
    pdf.set_text_color(80, 80, 80)
    pdf.cell(0, 7, f"Construct: {construct_name}",
             align="C", new_x="LMARGIN", new_y="NEXT")
    pdf.cell(0, 7,
             f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
             align="C", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(6)

    # ===== I. Construct Summary =====
    pdf.section_title("I. Construct Summary")

    summary_rows = [
        ["Insert gene", insert_name, f"{insert_len_bp} bp"],
        ["Vector backbone", vector_name, f"{vector_len_bp} bp"],
        ["Total construct", construct_name, f"{total_len} bp"],
    ]
    if selected_enzyme:
        summary_rows.append([
            "Restriction enzyme", selected_enzyme,
            f"Position: {insert_pos} bp",
        ])
    if protein_seq:
        detail = (f"{mw_kda:.1f} kDa / pI {pi_value:.2f}"
                  if mw_kda else "")
        summary_rows.append([
            "Encoded protein", f"{len(protein_seq)} AA", detail,
        ])
    if plddt is not None:
        confidence = ("High" if plddt > 70
                      else ("Medium" if plddt > 50 else "Low"))
        summary_rows.append([
            "ESMFold pLDDT", f"{plddt:.1f}",
            f"Confidence: {confidence}",
        ])

    pdf.add_table(["Item", "Value", "Detail"], summary_rows,
                  col_widths=[55, 70, 65])

    # Vector map image
    if vector_map_fig is not None:
        import tempfile
        import os
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            vector_map_fig.savefig(
                tmp.name, format="png", dpi=150, bbox_inches="tight")
            tmp_path = tmp.name
        try:
            pdf.image(tmp_path, x=40, w=130)
            pdf.ln(5)
        finally:
            os.unlink(tmp_path)

    # ===== II. PCR Protocol =====
    pdf.add_page()
    pdf.section_title(f"II. PCR Protocol ({proto['name']})")

    pdf.section_title("A. Reaction Mix", level=2)
    pcr_rows = [
        [proto["premix"], f"{proto['premix_vol']} uL", "1x"],
        ["Primer Forward",
         f"{proto['primer_pmol']} pmol", proto["primer_conc"]],
        ["Primer Reverse",
         f"{proto['primer_pmol']} pmol", proto["primer_conc"]],
        ["Template DNA", proto["template_max"], "-"],
        ["Sterile H2O", f"up to {proto['total_vol']} uL", "-"],
        ["Total", f"{proto['total_vol']} uL", ""],
    ]
    pdf.add_table(["Reagent", "Volume", "Final Conc."], pcr_rows,
                  col_widths=[80, 50, 60])

    pdf.note_box(
        "Reaction mix can be prepared at room temperature, "
        "but keep enzyme stocks on ice.", "note")

    # Tm / annealing time
    anneal_time = proto["anneal_time_high_tm"]
    tm_info = ""
    if primer_f:
        tm_f = calc_tm_simple(primer_f)
        tm_r = calc_tm_simple(primer_r) if primer_r else tm_f
        anneal_time = (proto["anneal_time_high_tm"]
                       if min(tm_f, tm_r) >= 55
                       else proto["anneal_time_low_tm"])
        tm_info = f"Estimated Tm: Fwd={tm_f:.0f} C / Rev={tm_r:.0f} C"

    pdf.section_title("B. Cycling Conditions", level=2)
    pdf.body_text(f"Target amplicon: {ext_times['insert_kb']} kb")
    if tm_info:
        pdf.body_text(tm_info)

    # Fast conditions
    pdf.body_text(
        f"(A) Standard (template < 200 ng / {proto['total_vol']} uL):")
    cycling_a = [
        [f"{proto['denature_temp']} C",
         f"{proto['denature_time']} sec", "Denaturation"],
        [f"{proto['anneal_temp']} C",
         f"{anneal_time} sec", "Annealing"],
        [f"{proto['extension_temp']} C",
         f"{ext_times['fast_sec']} sec",
         f"Extension ({proto['extension_rate_fast']} sec/kb)"],
    ]
    pdf.add_table(["Temp", "Time", "Step"], cycling_a,
                  col_widths=[40, 50, 100])
    pdf.body_text(
        f"Repeat for {proto['cycles_min']}"
        f" - {proto['cycles_max']} cycles.")

    # High-template conditions
    pdf.body_text(
        f"(B) High-template (> 200 ng / {proto['total_vol']} uL):")
    cycling_b = [
        [f"{proto['denature_temp']} C",
         f"{proto['denature_time']} sec", "Denaturation"],
        [f"{proto['anneal_temp']} C",
         f"{anneal_time} sec", "Annealing"],
        [f"{proto['extension_temp']} C",
         f"{ext_times['slow_sec']} sec",
         f"Extension ({proto['extension_rate_slow']} sec/kb)"],
    ]
    pdf.add_table(["Temp", "Time", "Step"], cycling_b,
                  col_widths=[40, 50, 100])

    pdf.note_box(
        "For PrimeSTAR Max, keep annealing time at 5 or 15 sec. "
        "Longer annealing may produce smears. "
        "If smears persist with 3-step, try 2-step PCR.", "warning")

    # ===== III. Restriction Enzyme Digestion =====
    if selected_enzyme:
        pdf.add_page()
        pdf.section_title("III. Restriction Enzyme Digestion")

        enz_info = ENZYME_CONDITIONS.get(selected_enzyme, {})
        if enz_info:
            digest_rows = [
                ["Enzyme", selected_enzyme],
                ["Recognition overhang",
                 enz_info.get("overhang", "-")],
                ["Recommended buffer",
                 enz_info.get("buffer", "-")],
                ["Reaction temp",
                 f"{enz_info.get('temp', 37)} C"],
                ["Incubation time", "1-2 hours"],
                ["Star activity risk",
                 enz_info.get("star_risk", "-")],
            ]
            pdf.add_table(["Parameter", "Value"], digest_rows,
                          col_widths=[80, 110])

        pdf.section_title("Digestion Reaction Mix", level=2)
        digest_mix = [
            ["DNA (vector or PCR product)", "1 ug", "-"],
            ["10x Restriction Buffer", "5 uL", "1x"],
            [f"{selected_enzyme}", "10-20 U", "-"],
            ["Sterile H2O", "up to 50 uL", "-"],
            ["Total", "50 uL", ""],
        ]
        pdf.add_table(["Reagent", "Amount", "Final Conc."], digest_mix,
                      col_widths=[80, 50, 60])

        if enz_info.get("overhang", "").startswith("3'"):
            pdf.note_box(
                f"{selected_enzyme} produces 3' overhangs. "
                "If using PrimeSTAR Max, perform protein removal "
                "(phenol/chloroform or PCR cleanup) BEFORE digestion "
                "to prevent 3'->5' exonuclease trimming of overhangs.",
                "warning")

        pdf.note_box(
            "After digestion, purify fragments using NucleoSpin "
            "Gel and PCR Clean-up or equivalent column purification. "
            "Gel-extract if multiple bands are present.", "tip")

    # ===== IV. Ligation =====
    pdf.add_page()
    pdf.section_title("IV. Ligation")

    is_blunt = (
        selected_enzyme
        and ENZYME_CONDITIONS.get(selected_enzyme, {})
            .get("overhang", "").startswith("Blunt")
    )
    lig_key = "blunt_end" if is_blunt else "sticky_end"
    lig = LIGATION_PROTOCOLS[lig_key]

    pdf.body_text(f"Protocol type: {lig['name']}")

    lig_rows = [
        ["10x T4 DNA Ligase Buffer", "2 uL", "1x"],
        ["Linearized vector", "50-100 ng", "-"],
        [f"Insert DNA ({lig['vector_insert_ratio']})",
         "Calculated", "-"],
        ["T4 DNA Ligase", "1 uL (200-400 U)", "-"],
        ["Sterile H2O", "up to 20 uL", "-"],
        ["Total", "20 uL", ""],
    ]
    pdf.add_table(["Reagent", "Amount", "Final Conc."], lig_rows,
                  col_widths=[80, 50, 60])

    # Insert amount calculation
    if insert_len_bp and vector_len_bp:
        pdf.section_title("Insert Amount Calculation", level=2)
        pdf.body_text(
            "Insert (ng) = [Vector (ng) x Insert size (kb) "
            "/ Vector size (kb)] x Molar ratio")
        ratio = 3
        insert_ng = (100 * insert_len_bp / vector_len_bp) * ratio
        pdf.body_text(
            f"Example (vector=100 ng, ratio={ratio}:1): "
            f"Insert = (100 x {insert_len_bp/1000:.2f} "
            f"/ {vector_len_bp/1000:.2f}) x {ratio} "
            f"= {insert_ng:.1f} ng")

    pdf.body_text(
        f"Incubate at {lig['temp']} C for {lig['time']}.")
    pdf.note_box(
        "Include a vector-only (no insert) control "
        "to assess background.", "tip")

    # ===== V. Transformation =====
    pdf.section_title("V. Transformation")
    pdf.body_text(
        "1. Thaw competent cells (e.g., DH5alpha, TOP10) "
        "on ice for 10 min.\n"
        "2. Add 2-5 uL of ligation mix to 50 uL competent cells. "
        "Mix gently.\n"
        "3. Incubate on ice for 30 min.\n"
        "4. Heat shock at 42 C for 45 sec "
        "(or per manufacturer protocol).\n"
        "5. Return to ice for 2 min.\n"
        "6. Add 450 uL SOC medium (pre-warmed to 37 C).\n"
        "7. Incubate at 37 C with shaking (200 rpm) for 1 hour.\n"
        "8. Plate 100-200 uL on selective agar plates.\n"
        "9. Incubate overnight at 37 C.")
    pdf.note_box(
        "Plate remaining cells by brief centrifugation "
        "and resuspension if low colony count is expected.", "tip")

    # ===== VI. Colony Screening =====
    pdf.add_page()
    pdf.section_title("VI. Colony Screening")
    pdf.body_text(
        "Pick 6-12 colonies with sterile tips into:\n"
        "  (a) 20 uL colony PCR reaction, AND\n"
        "  (b) numbered grid plate with selective medium.\n\n"
        "Colony PCR conditions: use vector-specific primers "
        "flanking the MCS.\n"
        f"Expected band: ~{insert_len_bp} bp (insert) "
        f"+ flanking vector region.")
    pdf.body_text(
        "Positive clones: miniprep and confirm by restriction "
        "digestion and/or sequencing.")

    # ===== VII. Predicted Protein Properties =====
    if protein_seq:
        pdf.section_title("VII. Predicted Protein Properties")
        prop_rows = [
            ["Amino acid length", f"{len(protein_seq)} AA"],
            ["Molecular weight",
             f"{mw_kda:.2f} kDa" if mw_kda else "-"],
            ["Isoelectric point (pI)",
             f"{pi_value:.2f}" if pi_value else "-"],
            ["ESMFold avg pLDDT",
             f"{plddt:.1f}" if plddt else "Not computed"],
        ]
        pdf.add_table(["Property", "Value"], prop_rows,
                      col_widths=[80, 110])

        if plddt is not None:
            if plddt > 70:
                pdf.note_box(
                    "High pLDDT (>70): Predicted structure is likely "
                    "reliable. Good candidate for expression.", "tip")
            elif plddt > 50:
                pdf.note_box(
                    "Medium pLDDT (50-70): Some regions may be "
                    "disordered. Consider truncation constructs.",
                    "note")
            else:
                pdf.note_box(
                    "Low pLDDT (<50): Protein may be largely "
                    "disordered or the prediction is unreliable. "
                    "Reconsider construct design.", "warning")

    # ===== VIII. Troubleshooting =====
    pdf.add_page()
    pdf.section_title("VIII. Troubleshooting")

    ts_rows = [
        ["No PCR product", "Ext. time too short",
         f"Set ext. to {ext_times['slow_sec']} sec or longer"],
        ["No PCR product", "Low template",
         "Increase to 35-40 cycles"],
        ["No PCR product", "Annealing issue",
         "Try 15 sec anneal / lower temp to 50-53 C"],
        ["Smear/extra bands", "Annealing too long",
         "Reduce to 5 sec"],
        ["Smear/extra bands", "Non-specific priming",
         "Raise anneal temp to 58-63 C / 2-step PCR"],
        ["Smear/extra bands", "Excess template",
         "Reduce template amount"],
        ["No colonies", "Ligation failure",
         "Check insert/vector ratio; use fresh ligase"],
        ["No colonies", "Dead competent cells",
         "Test with uncut plasmid control"],
        ["All empty colonies", "Incomplete digestion",
         "Re-digest vector; use phosphatase (CIP/SAP)"],
    ]
    pdf.add_table(["Problem", "Possible Cause", "Solution"], ts_rows,
                  col_widths=[45, 50, 95])

    return pdf.output()
