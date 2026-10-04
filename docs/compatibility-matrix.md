<!-- GENERATED FILE -- do not edit.
Regenerate with:  python3 oracle/gen_compat_matrix.py --out docs/compatibility-matrix.md
Source of truth: tests/golden/**/command.json (the recorded invocations)
-->

# Compatibility matrix

Generated from **7685 recorded invocations** by `oracle/gen_compat_matrix.py`. Each cell is `byte-identical files / recorded files` across the fixtures for that (subcommand, option variant). A dash means the pair is not covered by the recorded corpus, which is a gap in the matrix, not a pass.

| subcommand | B | broad | bw300 | call_summits | default | gsize_numeric | keepdup1 | keepdup_all | keepdup_auto | mfold_3_20 | mfold_bad_arity | nolambda | nomodel_extsize | nomodel_shift | q001 | q05 | scale_to_large | shift_only | slocal_500_llocal_2000 | spmr |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `callpeak` | 1256/1266 | 1254/1266 | 1256/1266 | 1247/1266 | 1256/1266 | 1256/1266 | 773/777 | 774/780 | 773/777 | 1256/1266 | 0/0 | 769/777 | 1256/1266 | 1252/1266 | 1256/1266 | 1256/1266 | 1250/1266 | 1252/1266 | 767/777 | 1256/1266 |

**Totals: 21415/21612 compared output files byte-identical** across 20 (subcommand, variant) pairs. 0 further recorded files were not produced by the replay and are excluded from the denominator.

> Not all recorded files are byte-identical yet; see `docs/upstream-findings.md` for the open findings.
