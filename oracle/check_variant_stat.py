#!/usr/bin/env python
"""Generate `crates/macs-callvar/tests/data/variant_stat.golden` from the oracle.

`MACS3/Signal/VariantStat.py` is where `callvar`'s allele likelihood lives. It is
compiled Cython, so the `.py` source is only a rough guide -- in particular the
`cython.float` parameter and the C `(int)` casts are decided by the compiler, not
by the text. This drives the *installed `.so`* with a few hundred cases spanning
the shapes the implementation branches on:

* both samples empty (everything is zero);
* one sample empty (control only, treatment only);
* a single read (`tn == 1`, which takes its own branch);
* every read top1 (`m == tn`);
* no read top1 (`m == 0`);
* and randomised bulk, including deliberately extreme base qualities.

Every numeric field is emitted with `%.17g` so the comparison is bit-exact rather
than tolerance-based: these functions are summed in a fixed order, and a `f64`
implementation that reassociates them would be a real behavioural difference even
though the peaks it feeds are identical.

    python3 oracle/check_variant_stat.py [--out PATH] [--cases N]
"""

import argparse
import datetime
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ENV_LOCK = os.path.join(HERE, "ENV.lock")


def load_env_lock():
    """Run the pinned oracle's Python so the comparison uses *its* NumPy/`.so`.

    The root is derived from `ENV.lock`'s `MACS3_PATH` the same way
    `oracle/check_predictd_randsample.sh` does -- `MACS3_PATH` points at
    `.../MACS3/__init__.py`, so the importable root is two levels up.
    """
    env = dict(os.environ)
    macs3_init = None
    with open(ENV_LOCK) as fh:
        for line in fh:
            if line.startswith("MACS3_PATH="):
                macs3_init = line.split("=", 1)[1].strip()
    if not macs3_init:
        raise SystemExit(f"{ENV_LOCK} has no MACS3_PATH")
    oracle_src = os.path.dirname(os.path.dirname(macs3_init))
    env["PYTHONPATH"] = oracle_src
    return env, oracle_src


def gen_cases(rng, n_bulk):
    """Yield (name, top1_T, top1_C, top2_T, top2_C) tuples."""
    yield ("empty_all", [], [], [], [])
    yield ("ctrl_only", [], [30, 31, 28], [], [25, 27])
    yield ("treat_only", [30, 33, 31], [], [22, 24], [])
    yield ("single_read", [37], [], [], [])
    yield ("single_read_ctrl", [], [37], [], [])
    yield ("all_top1", [30, 31, 29, 35, 33], [30, 32, 28], [], [])
    yield ("all_top2", [], [], [30, 31, 29, 35, 33], [30, 32, 28])
    yield ("m_is_one", [40], [38, 39], [20, 21, 22], [19, 23])
    yield ("max_qual", [93, 93, 93], [93, 93], [0, 0], [0])
    yield ("zero_qual", [0, 0], [0], [0, 0], [0])

    for i in range(n_bulk):
        mt = rng.integers(1, 40)
        mc = rng.integers(0, 40)
        t1 = rng.integers(0, 94, size=mt).tolist()
        c1 = rng.integers(0, 94, size=mc).tolist()
        t2 = rng.integers(0, 94, size=int(mt * rng.random())).tolist()
        c2 = rng.integers(0, 94, size=int(mc * rng.random())).tolist()
        yield (f"bulk{i:04d}", t1, c1, t2, c2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(
        os.path.dirname(HERE), "crates", "macs-callvar", "tests", "data",
        "variant_stat.golden"))
    ap.add_argument("--cases", type=int, default=400)
    ap.add_argument("--max-ar", type=float, default=0.99)
    args = ap.parse_args()

    import numpy as np

    rng = np.random.default_rng(20240922)
    cases = list(gen_cases(rng, args.cases))

    # Run inside the oracle interpreter.
    body = f'''
import json, sys
import numpy as np
from MACS3.Signal.VariantStat import (CalModel_Homo, CalModel_Heter_noAS,
                                      CalModel_Heter_AS, calculate_GQ,
                                      calculate_GQ_heterASsig)
cases = json.loads(sys.stdin.read())
out = []
for name, t1, c1, t2, c2 in cases:
    t1 = np.array(t1, dtype=np.int32)
    c1 = np.array(c1, dtype=np.int32)
    t2 = np.array(t2, dtype=np.int32)
    c2 = np.array(c2, dtype=np.int32)
    rec = {{"name": name}}
    try:
        lnL, bic = CalModel_Homo(t1, c1, t2, c2)
        rec["homo"] = [float(lnL), float(bic)]
    except Exception as e:
        rec["homo"] = ["ERR", str(e)]
    try:
        lnL, bic = CalModel_Heter_noAS(t1, c1, t2, c2)
        rec["heter_noas"] = [float(lnL), float(bic)]
    except Exception as e:
        rec["heter_noas"] = ["ERR", str(e)]
    try:
        lnL, bic = CalModel_Heter_AS(t1, c1, t2, c2, float({args.max_ar!r}))
        rec["heter_as"] = [float(lnL), float(bic)]
    except Exception as e:
        rec["heter_as"] = ["ERR", str(e)]
    for tag, a, b, c in (("gq12", -1.0, -2.0, -3.0),
                         ("gq_eq", -5.0, -5.0, -5.0),
                         ("gq_far", -1.0, -80.0, -160.0),
                         ("gq_near", -0.001, -0.002, -0.0005)):
        rec[tag] = int(calculate_GQ(a, b, c))
    for tag, a, b in (("assig12", -1.0, -2.0),
                      ("assig_eq", -5.0, -5.0),
                      ("assig_far", -1.0, -90.0),
                      ("assig_near", -0.001, -0.002)):
        rec[tag] = int(calculate_GQ_heterASsig(a, b))
    out.append(rec)
print(json.dumps(out))
'''
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as fh:
        fh.write(body)
        script = fh.name
    try:
        env, oracle_src = load_env_lock()
        proc = subprocess.run([sys.executable, script],
                              input=__import__("json").dumps(cases),
                              capture_output=True, text=True, env=env)
        if proc.returncode != 0:
            raise SystemExit(f"oracle run failed:\n{proc.stderr[-4000:]}")
        payload = __import__("json").loads(proc.stdout)
    finally:
        os.unlink(script)
    # Each row carries its own inputs first, then that row's oracle answers, so the
    # Rust side can replay the case without re-deriving it.
    records = [(ins, rec) for ins, rec in zip(cases, payload)]

    stamp = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
    lines = [
        "# macs3-rs L4 golden -- MACS3/Signal/VariantStat.py via the pinned .so",
        f"# oracle source: {oracle_src}",
        f"# cases: {len(records)}   generated: {stamp}",
        "# all floats are %.17g and must match bit-for-bit (fixed-order summation)",
        "# fields: name top1_T top1_C top2_T top2_C (comma-joined base qualities) "
        "homo_lnL homo_BIC heter_noas_lnL heter_noas_BIC heter_as_lnL heter_as_BIC "
        "gq12 gq_eq gq_far gq_near assig12 assig_eq assig_far assig_near",
    ]
    for ins, r in records:
        def fmt(v):
            return "ERR" if isinstance(v, str) else "%.17g" % v

        def arr(v):
            return ",".join(str(int(x)) for x in v)

        homo = r["homo"]
        hnoas = r["heter_noas"]
        has = r["heter_as"]
        lines.append("\t".join([
            r["name"],
            arr(ins[1]), arr(ins[2]), arr(ins[3]), arr(ins[4]),
            fmt(homo[0]), fmt(homo[1]),
            fmt(hnoas[0]), fmt(hnoas[1]), fmt(has[0]), fmt(has[1]),
            str(r["gq12"]), str(r["gq_eq"]), str(r["gq_far"]), str(r["gq_near"]),
            str(r["assig12"]), str(r["assig_eq"]), str(r["assig_far"]),
            str(r["assig_near"]),
        ]))

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    n_err = sum(1 for _, r in records if r["heter_noas"][0] == "ERR")
    print(f"wrote {args.out}: {len(records)} cases ({n_err} empty-treatment ERRs)")


if __name__ == "__main__":
    main()