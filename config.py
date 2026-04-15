"""Vector2Fold 設定・定数モジュール"""

from Bio.Restriction import (
    RestrictionBatch, EcoRI, BamHI, HindIII, NotI, XhoI,
    SpeI, PstI, NcoI, NdeI, SalI, KpnI, SacI, SmaI
)

# --- 主要制限酵素 ---
MAJOR_ENZYMES = RestrictionBatch([
    EcoRI, BamHI, HindIII, NotI, XhoI, SpeI,
    PstI, NcoI, NdeI, SalI, KpnI, SacI, SmaI
])

# --- 制限酵素の反応条件データベース ---
# NOTE: values are ASCII-only for fpdf2 Helvetica compatibility in PDF output
ENZYME_CONDITIONS = {
    "EcoRI":  {"buffer": "H / EcoRI buffer", "temp": 37, "overhang": "5' overhang (AATT)", "star_risk": "Medium"},
    "BamHI":  {"buffer": "K / BamHI buffer", "temp": 37, "overhang": "5' overhang (GATC)", "star_risk": "High"},
    "HindIII":{"buffer": "M / HindIII buffer","temp": 37, "overhang": "5' overhang (AGCT)", "star_risk": "Medium"},
    "NotI":   {"buffer": "H", "temp": 37, "overhang": "5' overhang (GGCC)", "star_risk": "Low"},
    "XhoI":   {"buffer": "H / D", "temp": 37, "overhang": "5' overhang (TCGA)", "star_risk": "Medium"},
    "SpeI":   {"buffer": "M", "temp": 37, "overhang": "5' overhang (CTAG)", "star_risk": "Low"},
    "PstI":   {"buffer": "H / M", "temp": 37, "overhang": "3' overhang (TGCA)", "star_risk": "Medium"},
    "NcoI":   {"buffer": "K / H", "temp": 37, "overhang": "5' overhang (CATG)", "star_risk": "Low"},
    "NdeI":   {"buffer": "K / M", "temp": 37, "overhang": "5' overhang (TA)", "star_risk": "Low"},
    "SalI":   {"buffer": "H", "temp": 37, "overhang": "5' overhang (TCGA)", "star_risk": "High"},
    "KpnI":   {"buffer": "L", "temp": 37, "overhang": "3' overhang (GTAC)", "star_risk": "Low"},
    "SacI":   {"buffer": "L", "temp": 37, "overhang": "3' overhang (AGCT)", "star_risk": "Low"},
    "SmaI":   {"buffer": "L / SmaI buffer", "temp": 25, "overhang": "Blunt end", "star_risk": "Low"},
}

# --- ESMFold API ---
ESMFOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
ESMFOLD_MAX_LENGTH = 1000
ESMFOLD_TIMEOUT = 120

# --- PCR protocol template parameters ---
# NOTE: values are ASCII-only for fpdf2 Helvetica compatibility in PDF output
PCR_PROTOCOLS = {
    "PrimeSTAR_Max": {
        "name": "PrimeSTAR Max DNA Polymerase",
        "product_code": "R045A",
        "premix": "PrimeSTAR Max Premix (2x)",
        "premix_vol": 25,
        "total_vol": 50,
        "primer_pmol": "10-15",
        "primer_conc": "0.2-0.3 uM",
        "template_max": "< 200 ng",
        "denature_temp": 98,
        "denature_time": 10,
        "anneal_temp": 55,
        "anneal_time_high_tm": 5,
        "anneal_time_low_tm": 15,
        "extension_temp": 72,
        "extension_rate_fast": 5,  # sec/kb
        "extension_rate_slow": 30,  # sec/kb for high template
        "cycles_min": 30,
        "cycles_max": 35,
    },
    "standard_Taq": {
        "name": "Taq DNA Polymerase (Standard)",
        "premix": "10x PCR Buffer",
        "total_vol": 50,
        "denature_temp": 94,
        "denature_time": 30,
        "anneal_temp": 55,
        "anneal_time_high_tm": 30,
        "anneal_time_low_tm": 30,
        "extension_temp": 72,
        "extension_rate_fast": 60,  # sec/kb
        "cycles_min": 25,
        "cycles_max": 35,
    }
}

# --- Ligation conditions ---
LIGATION_PROTOCOLS = {
    "sticky_end": {
        "name": "Sticky-end Ligation",
        "enzyme": "T4 DNA Ligase",
        "buffer": "T4 DNA Ligase Buffer (1x)",
        "temp": 16,
        "time": "30 min - overnight",
        "vector_insert_ratio": "1:3 - 1:5 (molar ratio)",
    },
    "blunt_end": {
        "name": "Blunt-end Ligation",
        "enzyme": "T4 DNA Ligase",
        "buffer": "T4 DNA Ligase Buffer (1x)",
        "temp": 16,
        "time": "overnight (16-18 hours)",
        "vector_insert_ratio": "1:5 - 1:10 (molar ratio)",
    }
}
