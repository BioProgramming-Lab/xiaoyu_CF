"""Smoke test for the self-contained xiaoyu_CF package.

Run from the repository root:

    python test/package_smoke_test.py

or from this directory:

    python package_smoke_test.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import xiaoyu_CF as flow  # noqa: E402
import cytoflow.utility as util  # noqa: E402


DATA_DIR = REPO_ROOT / "test" / "26Mar2025" / "HuaiDan-New Experiment-20250326-1722"
SUMMARY_FILE = REPO_ROOT / "test" / "package_smoke_summary.json"


def positive_fraction(experiment, channel: str, low: float, high: float = 1e10) -> float:
    gated = flow.assist_gate(experiment, channel=channel, low=low, high=high)
    return len(gated) / len(experiment) if len(experiment) else 0.0


def main() -> None:
    assert DATA_DIR.exists(), f"Missing test data directory: {DATA_DIR}"
    assert not list(REPO_ROOT.rglob("*.so")), "Compiled .so files are still present"
    assert sorted(util.scale._scale_mapping.keys()) == ["linear", "log"]

    expr = flow.xiaoyu_Expr(str(DATA_DIR))
    expr.add_conc_condition(1e-5, dilu_dir="down")

    cell_type = flow.auto_gate_FSC_SSC(expr.expr, keep=0.6)
    de_doublets = flow.auto_gate_FSC_A_H(cell_type, keep=0.7)

    required_channels = [
        "FSC 488/10-A",
        "SSC 488/10-A",
        "FSC 488/10-H",
        "Alexa 647-A",
        "mCherry-A",
        "TagBFP-A",
        "549/15-488 nm citrine-A",
    ]
    missing = [channel for channel in required_channels if channel not in de_doublets.data]
    assert not missing, f"Missing expected channels: {missing}"

    fractions = {}
    for char in ["A", "B", "C"]:
        mix = flow.subset_by_char(de_doublets, char)
        epcam = flow.assist_gate(mix, channel="mCherry-A", low=5e6, high=1e9)
        egfr = flow.assist_gate(mix, channel="TagBFP-A", low=4e6, high=1e9)
        coexpression = flow.assist_gate(
            mix,
            channel="549/15-488 nm citrine-A",
            low=4e6,
            high=1e9,
        )

        fractions[char] = {
            "CHO_EGFR": positive_fraction(egfr, "Alexa 647-A", low=3e8),
            "CHO_EpCAM": positive_fraction(epcam, "Alexa 647-A", low=3e8),
            "CHO_EGFR_EpCAM": positive_fraction(coexpression, "Alexa 647-A", low=3e8),
        }

    medians = (
        de_doublets.data.groupby("sample", observed=False)["Alexa 647-A"]
        .median()
        .to_dict()
    )

    summary = {
        "data_dir": str(DATA_DIR.relative_to(REPO_ROOT)),
        "scale_mapping": sorted(util.scale._scale_mapping.keys()),
        "fcs_files": len(expr.filelist),
        "events_imported": int(len(expr.expr)),
        "events_after_fsc_ssc_gate": int(len(cell_type)),
        "events_after_doublet_gate": int(len(de_doublets)),
        "positive_fraction_by_mix": fractions,
        "alexa_647_median_by_sample": {str(k): float(v) for k, v in medians.items()},
    }

    SUMMARY_FILE.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"Wrote {SUMMARY_FILE}")


if __name__ == "__main__":
    main()
