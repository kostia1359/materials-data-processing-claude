import numpy as np
import pandas as pd

from qc import assign as A
from qc import stats as S
from qc.config import load_config
from qc.kpis.registry import SHORTLIST

CFG = load_config()
# brief defaults (the deployed config raises them to remove baseline false alarms; see DECISIONS.md)
CFG["zones"] = {"investigate": 2.5, "reject": 4.0, "multi_count_z": 2.5, "multi_count_n": 3, "investigate_count_z": 2.0, "investigate_count_n": 2}


def _frame(rng):
    rows = []
    for b, shift in (("3", 0.0), ("1", 3.0), ("2", -3.0)):
        for i in range(5):
            r = {"sample_id": f"s{b}{i}", "batch": b}
            for j, k in enumerate(SHORTLIST):
                base = 10.0 + j
                r[k] = base + (shift if j < 3 else 0) + rng.normal(0, 0.3)
                r[f"{k}_se"] = 0.1
            rows.append(r)
    return pd.DataFrame(rows)


def test_robust_scale_floor():
    s = S.robust_stats(np.array([1.0, 1.0, 1.0]), np.array([0.0]), 0.05)
    assert s["scale"] == 0.05


def test_assignment_and_verdict():
    df = _frame(np.random.default_rng(0))
    st = S.baseline_stats(df, CFG)
    assert st["n"] == 5
    model = A.build_model(df, st, CFG)
    assert model["loio_accuracy_fixed_baseline"] == 1.0
    x = df[df["batch"] == "1"].iloc[0].to_dict()
    a = A.assign(x, [], model, st, CFG)
    assert a["assigned_batch"] == "1"
    assert abs(sum(a["probabilities"].values()) - 1) < 1e-9
    v = S.verdict(x, st, CFG, [1.0] * 5)
    assert v["decision"] in ("INVESTIGATE", "REJECT")
    assert v["p_floor"] == 1 / 6
    sigs = model["signatures"]
    assert set(sigs["1"]["signature"]) <= set(SHORTLIST[:3])


def test_conformal_floor_never_below():
    df = _frame(np.random.default_rng(1))
    st = S.baseline_stats(df, CFG)
    x = df.iloc[0].to_dict()
    x[SHORTLIST[0]] += 100
    v = S.verdict(x, st, CFG, [0.5, 0.6, 0.7])
    assert v["conformal_rank_p"] >= v["p_floor"] - 1e-12
