<!-- GENERATED FILE -- do not edit.
Regenerate with:  python3 oracle/gen_compat_matrix.py --out docs/compatibility-matrix.md
Source of truth: tests/golden/**/command.json (the recorded invocations)
-->

# Compatibility matrix

Generated from **7685 recorded invocations** by `oracle/gen_compat_matrix.py`. Each cell is `byte-identical files / recorded files` across the fixtures for that (subcommand, option variant). A dash means the pair is not covered by the recorded corpus, which is a gap in the matrix, not a pass.

| subcommand | B | broad | bw300 | call_summits | default | gsize_numeric | keepdup1 | keepdup_all | keepdup_auto | mfold_3_20 | mfold_bad_arity | nolambda | nomodel_extsize | nomodel_shift | q001 | q05 | scale_to_large | shift_only | slocal_500_llocal_2000 | spmr |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `callpeak` | 2110/2110 | 1266/1266 | 1266/1266 | 1266/1266 | 1266/1266 | 1266/1266 | 777/777 | 780/780 | 777/777 | 1266/1266 | 0/0 | 777/777 | 1266/1266 | 1266/1266 | 1266/1266 | 1266/1266 | 1266/1266 | 1266/1266 | 777/777 | 1266/1266 |

**Totals: 22456/22456 compared output files byte-identical** across 20 (subcommand, variant) pairs. 0 recorded files were not produced by the replay and remain in the denominator. Exit-status mismatches: 0.
