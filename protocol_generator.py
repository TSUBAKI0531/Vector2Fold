"""
Vector2Fold Protocol Auto-Generation Module
Generates bench protocol PDFs for cloning experiments.
"""

import math
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
    "\u301C": "-",   # wave dash
    "\uFF5E": "-",   # fullwidth tilde
    "\u00D7": "x",   # multiplication sign
    "\u00B5": "u",   # micro sign
    "\u03BC": "u",   # Greek mu
    "\u2212": "-",   # minus sign
    "\u2013": "-",   # en dash
    "\u2014": "-",   # em dash
    "\u2018": "'",   # left single quote
    "\u2019": "'",   # right single quote
    "\u201C": '"',   # left double quote
    "\u201D": '"',   # right double quote
    "\u2026": "...", # ellipsis
    "\u00B0": " ",   # degree sign
    "\u2192": "->",  # right arrow
    "\u2190": "<-",  # left arrow
    "\u2265": ">=",
    "\u2264": "<=",
}


def _s(text) -> str:
    """Sanitize text for fpdf2 Helvetica (Latin-1 safe)."""
    t = str(text) if not isinstance(text, str) else text
    for char, repl in _UNICODE_REPLACEMENTS.items():
        t = t.replace(char, repl)
    return "".join(c if ord(c) <= 0x00FF else "?" for c in t)


# =========================================================================
# Custom PDF class (NO cell/multi_cell overrides)
# =========================================================================

class ProtocolPDF(FPDF):

    def __init__(self):
        super().__init__()
        self.set_auto_page_break(auto=True, margin=20)

    def header(self):
        self.set_font("Helvetica", "B", 10)
        self.set_text_color(100, 100, 100)
        self.cell(0, 6, _s("Vector2Fold - Auto-generated Bench Protocol"),
                  align="R", new_x="LMARGIN", new_y="NEXT")
        self.set_draw_color(0, 102, 204)
        self.set_line_width(0.5)
        self.line(10, self.get_y(), 200, self.get_y())
        self.ln(4)

    def footer(self):
        self.set_y(-15)
        self.set_font("Helvetica", "I", 8)
        self.set_text_color(128, 128, 128)
        self.cell(0, 10, _s(f"Page {self.page_no()}/{{nb}}"), align="C")

    def section_title(self, title, level=1):
        if level == 1:
            self.set_font("Helvetica", "B", 14)
            self.set_text_color(0, 51, 102)
            self.set_fill_color(230, 240, 250)
            self.cell(0, 10, _s(title), fill=True,
                      new_x="LMARGIN", new_y="NEXT")
        elif level == 2:
            self.set_font("Helvetica", "B", 11)
            self.set_text_color(0, 80, 140)
            self.cell(0, 8, _s(title),
                      new_x="LMARGIN", new_y="NEXT")
        self.ln(2)

    def body_text(self, text):
        self.set_font("Helvetica", "", 10)
        self.set_text_color(30, 30, 30)
        self.multi_cell(0, 5, _s(text))
        self.ln(2)

    def note_box(self, text, box_type="note"):
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
        self.cell(0, 6, _s(f"  {labels.get(box_type, '[NOTE]')}"),
                  fill=True, new_x="LMARGIN", new_y="NEXT")
        self.set_font("Helvetica", "", 9)
        self.set_text_color(60, 60, 60)
        self.multi_cell(0, 5, _s(f"  {text}"), fill=True)
        self.ln(3)

    def add_table(self, headers, rows, col_widths=None):
        if col_widths is None:
            col_widths = [190 / len(headers)] * len(headers)

        self.set_font("Helvetica", "B", 9)
        self.set_fill_color(0, 80, 140)
        self.set_text_color(255, 255, 255)
        for i, h in enumerate(headers):
            self.cell(col_widths[i], 7, _s(h),
                      border=1, fill=True, align="C")
        self.ln()

        self.set_font("Helvetica", "", 9)
        self.set_text_color(30, 30, 30)
        for row_idx, row in enumerate(rows):
            if row_idx % 2 == 0:
                self.set_fill_color(245, 248, 252)
            else:
                self.set_fill_color(255, 255, 255)
            for i, val in enumerate(row):
                self.cell(col_widths[i], 6, _s(val),
                          border=1, fill=True, align="C")
            self.ln()
        self.ln(3)


# =========================================================================
# Helpers
# =========================================================================

def calc_tm_simple(primer_seq):
    p = primer_seq.upper().strip()
    return (2 * (p.count("A") + p.count("T"))
            + 4 * (p.count("G") + p.count("C")) - 5)


def calc_pcr_extension_time(insert_len_bp, protocol_key="PrimeSTAR_Max"):
    proto = PCR_PROTOCOLS[protocol_key]
    kb = insert_len_bp / 1000
    return {
        "fast_sec": max(math.ceil(kb * proto["extension_rate_fast"]), 5),
        "slow_sec": max(math.ceil(kb * proto.get("extension_rate_slow",
                                                  proto["extension_rate_fast"] * 6)), 10),
        "insert_kb": round(kb, 2),
    }


# =========================================================================
# Main PDF generation
# =========================================================================

def generate_protocol_pdf(
    construct_name,
    insert_name,
    vector_name,
    insert_len_bp,
    vector_len_bp,
    selected_enzyme=None,
    insert_pos=None,
    protein_seq=None,
    mw_kda=None,
    pi_value=None,
    plddt=None,
    primer_f=None,
    primer_r=None,
    protocol_key="PrimeSTAR_Max",
    vector_map_fig=None,
):
    pdf = ProtocolPDF()
    pdf.alias_nb_pages()
    pdf.add_page()

    proto = PCR_PROTOCOLS[protocol_key]
    ext = calc_pcr_extension_time(insert_len_bp, protocol_key)
    total_len = insert_len_bp + vector_len_bp

    # Title
    pdf.set_font("Helvetica", "B", 18)
    pdf.set_text_color(0, 51, 102)
    pdf.cell(0, 12, _s("Bench Protocol"),
             align="C", new_x="LMARGIN", new_y="NEXT")
    pdf.set_font("Helvetica", "", 11)
    pdf.set_text_color(80, 80, 80)
    pdf.cell(0, 7, _s(f"Construct: {construct_name}"),
             align="C", new_x="LMARGIN", new_y="NEXT")
    pdf.cell(0, 7,
             _s(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}"),
             align="C", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(6)

    # I. Construct Summary
    pdf.section_title("I. Construct Summary")
    rows = [
        ["Insert gene", insert_name, f"{insert_len_bp} bp"],
        ["Vector backbone", vector_name, f"{vector_len_bp} bp"],
        ["Total construct", construct_name, f"{total_len} bp"],
    ]
    if selected_enzyme:
        rows.append(["Restriction enzyme", selected_enzyme,
                      f"Position: {insert_pos} bp"])
    if protein_seq:
        d = f"{mw_kda:.1f} kDa / pI {pi_value:.2f}" if mw_kda else ""
        rows.append(["Encoded protein", f"{len(protein_seq)} AA", d])
    if plddt is not None:
        c = "High" if plddt > 70 else ("Medium" if plddt > 50 else "Low")
        rows.append(["ESMFold pLDDT", f"{plddt:.1f}", f"Confidence: {c}"])
    pdf.add_table(["Item", "Value", "Detail"], rows,
                  col_widths=[55, 70, 65])

    # Vector map
    if vector_map_fig is not None:
        import tempfile, os
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            vector_map_fig.savefig(tmp.name, format="png", dpi=150,
                                  bbox_inches="tight")
            tmp_path = tmp.name
        try:
            pdf.image(tmp_path, x=40, w=130)
            pdf.ln(5)
        finally:
            os.unlink(tmp_path)

    # II. PCR Protocol
    pdf.add_page()
    pdf.section_title(f"II. PCR Protocol ({proto['name']})")
    pdf.section_title("A. Reaction Mix", level=2)

    pcr_rows = [
        [proto.get("premix", "Enzyme + Buffer"),
         f"{proto.get('premix_vol', '25')} uL", "1x"],
        ["Primer Forward",
         f"{proto.get('primer_pmol', '10-15')} pmol",
         proto.get("primer_conc", "0.2-0.3 uM")],
        ["Primer Reverse",
         f"{proto.get('primer_pmol', '10-15')} pmol",
         proto.get("primer_conc", "0.2-0.3 uM")],
        ["Template DNA", proto.get("template_max", "< 200 ng"), "-"],
        ["Sterile H2O", f"up to {proto['total_vol']} uL", "-"],
        ["Total", f"{proto['total_vol']} uL", ""],
    ]
    pdf.add_table(["Reagent", "Volume", "Final Conc."], pcr_rows,
                  col_widths=[80, 50, 60])
    pdf.note_box("Reaction mix can be prepared at room temperature, "
                 "but keep enzyme stocks on ice.", "note")

    # Annealing time
    anneal_time = proto["anneal_time_high_tm"]
    tm_info = ""
    if primer_f:
        tm_f = calc_tm_simple(primer_f)
        tm_r = calc_tm_simple(primer_r) if primer_r else tm_f
        anneal_time = (proto["anneal_time_high_tm"] if min(tm_f, tm_r) >= 55
                       else proto["anneal_time_low_tm"])
        tm_info = f"Estimated Tm: Fwd={tm_f:.0f} C / Rev={tm_r:.0f} C"

    pdf.section_title("B. Cycling Conditions", level=2)
    pdf.body_text(f"Target amplicon: {ext['insert_kb']} kb")
    if tm_info:
        pdf.body_text(tm_info)

    pdf.body_text(f"(A) Standard (template < 200 ng / {proto['total_vol']} uL):")
    pdf.add_table(["Temp", "Time", "Step"], [
        [f"{proto['denature_temp']} C", f"{proto['denature_time']} sec",
         "Denaturation"],
        [f"{proto['anneal_temp']} C", f"{anneal_time} sec", "Annealing"],
        [f"{proto['extension_temp']} C", f"{ext['fast_sec']} sec",
         f"Extension ({proto['extension_rate_fast']} sec/kb)"],
    ], col_widths=[40, 50, 100])
    pdf.body_text(f"Repeat for {proto['cycles_min']}"
                  f" - {proto['cycles_max']} cycles.")

    pdf.body_text(f"(B) High-template (> 200 ng / {proto['total_vol']} uL):")
    pdf.add_table(["Temp", "Time", "Step"], [
        [f"{proto['denature_temp']} C", f"{proto['denature_time']} sec",
         "Denaturation"],
        [f"{proto['anneal_temp']} C", f"{anneal_time} sec", "Annealing"],
        [f"{proto['extension_temp']} C", f"{ext['slow_sec']} sec",
         f"Extension ({proto.get('extension_rate_slow', 30)} sec/kb)"],
    ], col_widths=[40, 50, 100])
    pdf.note_box("For PrimeSTAR Max, keep annealing time at 5 or 15 sec. "
                 "Longer annealing may produce smears. "
                 "If smears persist with 3-step, try 2-step PCR.", "warning")

    # III. Restriction Enzyme Digestion
    if selected_enzyme:
        pdf.add_page()
        pdf.section_title("III. Restriction Enzyme Digestion")
        enz = ENZYME_CONDITIONS.get(selected_enzyme, {})
        if enz:
            pdf.add_table(["Parameter", "Value"], [
                ["Enzyme", selected_enzyme],
                ["Recognition overhang", enz.get("overhang", "-")],
                ["Recommended buffer", enz.get("buffer", "-")],
                ["Reaction temp", f"{enz.get('temp', 37)} C"],
                ["Incubation time", "1-2 hours"],
                ["Star activity risk", enz.get("star_risk", "-")],
            ], col_widths=[80, 110])

        pdf.section_title("Digestion Reaction Mix", level=2)
        pdf.add_table(["Reagent", "Amount", "Final Conc."], [
            ["DNA (vector or PCR product)", "1 ug", "-"],
            ["10x Restriction Buffer", "5 uL", "1x"],
            [selected_enzyme, "10-20 U", "-"],
            ["Sterile H2O", "up to 50 uL", "-"],
            ["Total", "50 uL", ""],
        ], col_widths=[80, 50, 60])

        if enz.get("overhang", "").startswith("3'"):
            pdf.note_box(
                f"{selected_enzyme} produces 3' overhangs. "
                "If using PrimeSTAR Max, perform protein removal "
                "(phenol/chloroform or PCR cleanup) BEFORE digestion "
                "to prevent 3'->5' exonuclease trimming.", "warning")
        pdf.note_box(
            "After digestion, purify fragments using NucleoSpin "
            "Gel and PCR Clean-up or equivalent. "
            "Gel-extract if multiple bands are present.", "tip")

    # IV. Ligation
    pdf.add_page()
    pdf.section_title("IV. Ligation")
    is_blunt = (selected_enzyme and
                ENZYME_CONDITIONS.get(selected_enzyme, {})
                .get("overhang", "").startswith("Blunt"))
    lig = LIGATION_PROTOCOLS["blunt_end" if is_blunt else "sticky_end"]
    pdf.body_text(f"Protocol type: {lig['name']}")
    pdf.add_table(["Reagent", "Amount", "Final Conc."], [
        ["10x T4 DNA Ligase Buffer", "2 uL", "1x"],
        ["Linearized vector", "50-100 ng", "-"],
        [f"Insert DNA ({lig['vector_insert_ratio']})", "Calculated", "-"],
        ["T4 DNA Ligase", "1 uL (200-400 U)", "-"],
        ["Sterile H2O", "up to 20 uL", "-"],
        ["Total", "20 uL", ""],
    ], col_widths=[80, 50, 60])

    if insert_len_bp and vector_len_bp:
        pdf.section_title("Insert Amount Calculation", level=2)
        pdf.body_text(
            "Insert (ng) = [Vector (ng) x Insert size (kb) "
            "/ Vector size (kb)] x Molar ratio")
        ratio = 3
        ins_ng = (100 * insert_len_bp / vector_len_bp) * ratio
        pdf.body_text(
            f"Example (vector=100 ng, ratio={ratio}:1): "
            f"Insert = (100 x {insert_len_bp/1000:.2f} "
            f"/ {vector_len_bp/1000:.2f}) x {ratio} = {ins_ng:.1f} ng")

    pdf.body_text(f"Incubate at {lig['temp']} C for {lig['time']}.")
    pdf.note_box("Include a vector-only (no insert) control "
                 "to assess background.", "tip")

    # V. Transformation
    pdf.section_title("V. Transformation")
    pdf.body_text(
        "1. Thaw competent cells (DH5alpha, TOP10) on ice 10 min.\n"
        "2. Add 2-5 uL ligation mix to 50 uL cells. Mix gently.\n"
        "3. Incubate on ice for 30 min.\n"
        "4. Heat shock at 42 C for 45 sec.\n"
        "5. Return to ice for 2 min.\n"
        "6. Add 450 uL SOC medium (pre-warmed 37 C).\n"
        "7. Incubate 37 C, shaking 200 rpm, 1 hour.\n"
        "8. Plate 100-200 uL on selective plates.\n"
        "9. Incubate overnight at 37 C.")
    pdf.note_box("Plate remaining cells by centrifugation "
                 "if low colony count expected.", "tip")

    # VI. Colony Screening
    pdf.add_page()
    pdf.section_title("VI. Colony Screening")
    pdf.body_text(
        "Pick 6-12 colonies with sterile tips into:\n"
        "  (a) 20 uL colony PCR reaction, AND\n"
        "  (b) numbered grid plate with selective medium.\n\n"
        "Colony PCR: use vector-specific primers flanking MCS.\n"
        f"Expected band: ~{insert_len_bp} bp + flanking region.")
    pdf.body_text(
        "Positive clones: miniprep and confirm by digestion "
        "and/or sequencing.")

    # VII. Protein Properties
    if protein_seq:
        pdf.section_title("VII. Predicted Protein Properties")
        pdf.add_table(["Property", "Value"], [
            ["Amino acid length", f"{len(protein_seq)} AA"],
            ["Molecular weight",
             f"{mw_kda:.2f} kDa" if mw_kda else "-"],
            ["Isoelectric point (pI)",
             f"{pi_value:.2f}" if pi_value else "-"],
            ["ESMFold avg pLDDT",
             f"{plddt:.1f}" if plddt else "Not computed"],
        ], col_widths=[80, 110])

        if plddt is not None:
            if plddt > 70:
                pdf.note_box(
                    "High pLDDT (>70): Structure likely reliable. "
                    "Good candidate for expression.", "tip")
            elif plddt > 50:
                pdf.note_box(
                    "Medium pLDDT (50-70): Some disordered regions. "
                    "Consider truncation constructs.", "note")
            else:
                pdf.note_box(
                    "Low pLDDT (<50): Largely disordered or unreliable. "
                    "Reconsider construct design.", "warning")

    # VIII. Troubleshooting
    pdf.add_page()
    pdf.section_title("VIII. Troubleshooting")
    pdf.add_table(["Problem", "Possible Cause", "Solution"], [
        ["No PCR product", "Ext. time too short",
         f"Set ext. to {ext['slow_sec']} sec or longer"],
        ["No PCR product", "Low template",
         "Increase to 35-40 cycles"],
        ["No PCR product", "Annealing issue",
         "Try 15 sec / lower temp to 50-53 C"],
        ["Smear/extra bands", "Annealing too long",
         "Reduce to 5 sec"],
        ["Smear/extra bands", "Non-specific priming",
         "Raise temp to 58-63 C / 2-step PCR"],
        ["Smear/extra bands", "Excess template",
         "Reduce template amount"],
        ["No colonies", "Ligation failure",
         "Check ratio; use fresh ligase"],
        ["No colonies", "Dead comp. cells",
         "Test with uncut plasmid control"],
        ["All empty colonies", "Incomplete digestion",
         "Re-digest; use phosphatase (CIP/SAP)"],
    ], col_widths=[45, 50, 95])

    return bytes(pdf.output())
