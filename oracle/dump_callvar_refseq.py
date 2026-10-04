#!/usr/bin/env python
"""Dump upstream's own `RACollection` consensus per peak, for the callvar differential.

`RACollection` builds the peak's reference sequence *from the reads* -- there is no
genome FASTA -- and `__fill_refseq` fills it with a Python `bytearray` slice
assignment whose length differs from the slice it replaces. That resize is the thing
this port has to reproduce, and reproducing it correctly is hard to eyeball from a VCF
diff: a single shifted byte changes every reference base downstream.

So this dumps the *inputs and output* of that construction, per peak, straight from the
pinned oracle:

    chrom left right RAs_left RAs_right count_T count_C len(refseq) hex(refseq)

and `macs-callvar` produces the identical record for the same peak. The checker then
compares them field by field, which localises a divergence to a specific peak instead
of to "the VCF".

    python3 oracle/dump_callvar_refseq.py [--out PATH]
"""

import argparse
import datetime
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ENV_LOCK = os.path.join(HERE, "ENV.lock")


def load_env_lock():
    env = dict(os.environ)
    macs3_init = None
    with open(ENV_LOCK) as fh:
        for line in fh:
            if line.startswith("MACS3_PATH="):
                macs3_init = line.split("=", 1)[1].strip()
    if not macs3_init:
        raise SystemExit(f"{ENV_LOCK} has no MACS3_PATH")
    src = os.path.dirname(os.path.dirname(macs3_init))
    env["PYTHONPATH"] = src
    return env, src


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--peaks", default="/tmp/cv/callvar_testing.narrowPeak")
    ap.add_argument("--tbam", default="/tmp/cv/CTCF_PE_ChIP_chr22_50k.bam")
    ap.add_argument("--cbam", default="/tmp/cv/CTCF_PE_CTRL_chr22_50k.bam")
    ap.add_argument("--out", default=os.path.join(
        ROOT, "crates", "macs-callvar", "tests", "data",
        "callvar_peak_refseq.golden"))
    args = ap.parse_args()

    body = r'''
import json, sys
from MACS3.IO.BAM import BAMaccessor
from MACS3.Signal.RACollection import RACollection

peaks, tbam_p, cbam_p, maxdup = sys.argv[1:5]
maxdup = int(maxdup)
tbam = BAMaccessor(tbam_p)
cbam = BAMaccessor(cbam_p)

out = []
with open(peaks) as fh:
    for line in fh:
        f = line.rstrip().split()
        if len(f) < 3:
            continue
        chrom, start, end = f[0].encode(), int(f[1]), int(f[2])
        try:
            if cbam_p:
                coll = RACollection(chrom, {"start": start, "end": end},
                                     tbam.get_reads_in_region(chrom, start, end, maxDuplicate=maxdup),
                                     cbam.get_reads_in_region(chrom, start, end, maxDuplicate=maxdup))
            else:
                coll = RACollection(chrom, {"start": start, "end": end},
                                     tbam.get_reads_in_region(chrom, start, end, maxDuplicate=maxdup))
        except Exception:
            # "No reads found in this peak. Skipped"
            continue
        seq = bytes(coll["peak_refseq"])
        out.append({
            "chrom": chrom.decode(),
            "left": coll["left"], "right": coll["right"],
            "ras_left": coll["RAs_left"], "ras_right": coll["RAs_right"],
            "count_t": coll["count_T"], "count_c": coll["count_C"],
            "ext_len": len(bytes(coll["peak_refseq_ext"])),
            "seq": seq.hex(),
            "ext": bytes(coll["peak_refseq_ext"]).hex(),
        })
print(json.dumps(out))
'''

    env, src = load_env_lock()
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as fh:
        fh.write(body)
        script = fh.name
    try:
        proc = subprocess.run([sys.executable, script, args.peaks, args.tbam,
                               args.cbam, "1"],
                              capture_output=True, text=True, env=env, cwd=src)
        if proc.returncode != 0:
            raise SystemExit(f"oracle run failed:\n{proc.stderr[-4000:]}")
        import json
        recs = json.loads(proc.stdout)
    finally:
        os.unlink(script)

    stamp = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
    lines = [
        "# macs3-rs L4 golden -- upstream's own RACollection peak consensus.",
        f"# oracle source: {src}",
        f"# peaks: {len(recs)} (read-covered ones only)   generated: {stamp}",
        "# fields: chrom left right RAs_left RAs_right count_T count_C ext_len refseq_hex ext_hex",
        "# refseq_hex is peak_refseq, i.e. the slice [left-start : right-start] of",
        "# peak_refseq_ext AFTER the bytearray resizes the fill performed.",
    ]
    for r in recs:
        lines.append("\t".join(str(r[k]) for k in
                    ("chrom", "left", "right", "ras_left", "ras_right",
                     "count_t", "count_c", "ext_len", "seq", "ext")))
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"wrote {args.out}: {len(recs)} peaks")


if __name__ == "__main__":
    main()