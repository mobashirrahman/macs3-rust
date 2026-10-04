# fermi-lite, vendored

This is MACS3's own fermi-lite submodule, copied **verbatim** from the pinned oracle at
`oracle/ENV.lock`'s `MACS3_COMMIT` (`c5443190e3edfeb301cc94cc450e2b2c026a223`):

- upstream path: `MACS3/fermi-lite/`
- submodule commit recorded by the oracle's `.gitmodules`: `c1e71a7aec22fdbbf2f913a4bd1bcc352f7b99ac`
- fermi-lite reports its own version as `FML_VERSION "r53"`
- licence: MIT, Copyright (c) 2016 Broad Institute — see `LICENSE.txt`

## Why it is vendored rather than linked from the oracle tree

`PORTING_PLAN.md` permits exactly this arrangement: *"callvar may initially bridge to
upstream fermi-lite via thin C FFI; a pure-Rust assembler is the last item, after every
other gate is green."*

Vendoring means the build does not require the oracle checkout to be present, and — more
importantly — that the assembler is **the same code** upstream runs, so `--fermi auto`
can be made byte-identical rather than approximately right. A reimplementation of a de
novo assembler would be a research project with its own failure modes; the declared
deviation budget for this project is already spent on HMMRATAC self-training.

## What is and is not here

`example.c` is **not** vendored: it defines `main`, and only
`rle.c`, `rope.c`, `unitig.c` and `misc.c` are compiled (the same set fermi-lite's own
`Makefile` builds for the library).

Nothing in this directory has been edited. The only Rust-side adaptation is the build
script (`crates/macs-callvar/build.rs`) and the FFI declarations
(`crates/macs-callvar/src/fermi.rs`).

## Rebuilding / replacing

To move to a different fermi-lite revision, replace this directory wholesale and update
the commit above. Do not hand-edit the C: any divergence here is divergence from the
variant-calling oracle.