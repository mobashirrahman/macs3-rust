//! `-p/--pvalue` must threshold the **p-score** track, not the q-score track.
//!
//! `PeakDetect` dispatches on which cutoff the user gave: `call_peaks(['p'], [log_pvalue])`
//! when `-p` is set, `call_peaks(['q'], [log_qvalue])` otherwise (`PeakDetect.py:255-289`
//! narrow, `:372-406` without a control). One `score_array_s` is built per call from
//! that symbol (`CallPeakUnit.py:1184-1197`, and `:1983-1996` for both broad levels) and
//! thresholded by `apply_multiple_cutoffs(score_array_s, cutoff_s)`. This port always
//! compared the **q**-score against the p-value cutoff, so `-p 0.01` silently produced
//! the `-q 0.01` answer: on the 5 M-read fixture upstream calls 69,159 narrow peaks
//! under `-p 0.01` where this port emitted 35,683, and under `--broad -p 0.01` upstream
//! emits 1,132,451 rows against our 40,754.
//!
//! The fixture below is the smallest input that separates the two paths: with it the
//! reference calls three peaks under `-p 0.01` but only two under `-q 0.01`, so the
//! expected output cannot be produced by thresholding the wrong track. Every expected
//! string is a verbatim copy of the pinned oracle's output for this fixture
//! (`OPENBLAS_CORETYPE=Haswell .oracle/venv/bin/macs3 ...`); nothing here needs Python
//! or the oracle at run time.

const RUST: &str = env!("CARGO_BIN_EXE_macs3-rs");

/// 135 single-end reads in three clusters on a 20 kb "genome".
const TREAT: &str = concat!(
    "chr1	2000	2030	.	.	+
",
    "chr1	2000	2030	.	.	+
",
    "chr1	2000	2030	.	.	+
",
    "chr1	2001	2031	.	.	+
",
    "chr1	2003	2033	.	.	+
",
    "chr1	2004	2034	.	.	+
",
    "chr1	2006	2036	.	.	+
",
    "chr1	2007	2037	.	.	+
",
    "chr1	2008	2038	.	.	+
",
    "chr1	2008	2038	.	.	+
",
    "chr1	2008	2038	.	.	+
",
    "chr1	2010	2040	.	.	+
",
    "chr1	2010	2040	.	.	+
",
    "chr1	2010	2040	.	.	+
",
    "chr1	2010	2040	.	.	+
",
    "chr1	2011	2041	.	.	+
",
    "chr1	2011	2041	.	.	+
",
    "chr1	2012	2042	.	.	+
",
    "chr1	2013	2043	.	.	+
",
    "chr1	2013	2043	.	.	+
",
    "chr1	2013	2043	.	.	+
",
    "chr1	2013	2043	.	.	+
",
    "chr1	2015	2045	.	.	+
",
    "chr1	2015	2045	.	.	+
",
    "chr1	2015	2045	.	.	+
",
    "chr1	2016	2046	.	.	+
",
    "chr1	2017	2047	.	.	+
",
    "chr1	2018	2048	.	.	+
",
    "chr1	2020	2050	.	.	+
",
    "chr1	2022	2052	.	.	+
",
    "chr1	2023	2053	.	.	+
",
    "chr1	2024	2054	.	.	+
",
    "chr1	2024	2054	.	.	+
",
    "chr1	2026	2056	.	.	+
",
    "chr1	2028	2058	.	.	+
",
    "chr1	2029	2059	.	.	+
",
    "chr1	2030	2060	.	.	+
",
    "chr1	2033	2063	.	.	+
",
    "chr1	2034	2064	.	.	+
",
    "chr1	2034	2064	.	.	+
",
    "chr1	2036	2066	.	.	+
",
    "chr1	2039	2069	.	.	+
",
    "chr1	2039	2069	.	.	+
",
    "chr1	2039	2069	.	.	+
",
    "chr1	2040	2070	.	.	+
",
    "chr1	2600	2630	.	.	+
",
    "chr1	2601	2631	.	.	+
",
    "chr1	2601	2631	.	.	+
",
    "chr1	2603	2633	.	.	+
",
    "chr1	2604	2634	.	.	+
",
    "chr1	2604	2634	.	.	+
",
    "chr1	2609	2639	.	.	+
",
    "chr1	2610	2640	.	.	+
",
    "chr1	2611	2641	.	.	+
",
    "chr1	2611	2641	.	.	+
",
    "chr1	2612	2642	.	.	+
",
    "chr1	2616	2646	.	.	+
",
    "chr1	2616	2646	.	.	+
",
    "chr1	2619	2649	.	.	+
",
    "chr1	2619	2649	.	.	+
",
    "chr1	2619	2649	.	.	+
",
    "chr1	2619	2649	.	.	+
",
    "chr1	2620	2650	.	.	+
",
    "chr1	2621	2651	.	.	+
",
    "chr1	2621	2651	.	.	+
",
    "chr1	2622	2652	.	.	+
",
    "chr1	2623	2653	.	.	+
",
    "chr1	2624	2654	.	.	+
",
    "chr1	2626	2656	.	.	+
",
    "chr1	2630	2660	.	.	+
",
    "chr1	2630	2660	.	.	+
",
    "chr1	2630	2660	.	.	+
",
    "chr1	2637	2667	.	.	+
",
    "chr1	2638	2668	.	.	+
",
    "chr1	2638	2668	.	.	+
",
    "chr1	4200	4230	.	.	+
",
    "chr1	4201	4231	.	.	+
",
    "chr1	4202	4232	.	.	+
",
    "chr1	4202	4232	.	.	+
",
    "chr1	4205	4235	.	.	+
",
    "chr1	4205	4235	.	.	+
",
    "chr1	4206	4236	.	.	+
",
    "chr1	4207	4237	.	.	+
",
    "chr1	4209	4239	.	.	+
",
    "chr1	4211	4241	.	.	+
",
    "chr1	4212	4242	.	.	+
",
    "chr1	4213	4243	.	.	+
",
    "chr1	4215	4245	.	.	+
",
    "chr1	4216	4246	.	.	+
",
    "chr1	4217	4247	.	.	+
",
    "chr1	4218	4248	.	.	+
",
    "chr1	4221	4251	.	.	+
",
    "chr1	4221	4251	.	.	+
",
    "chr1	4222	4252	.	.	+
",
    "chr1	4222	4252	.	.	+
",
    "chr1	4222	4252	.	.	+
",
    "chr1	4223	4253	.	.	+
",
    "chr1	4223	4253	.	.	+
",
    "chr1	4223	4253	.	.	+
",
    "chr1	4224	4254	.	.	+
",
    "chr1	4225	4255	.	.	+
",
    "chr1	4226	4256	.	.	+
",
    "chr1	4227	4257	.	.	+
",
    "chr1	4228	4258	.	.	+
",
    "chr1	4229	4259	.	.	+
",
    "chr1	4229	4259	.	.	+
",
    "chr1	4232	4262	.	.	+
",
    "chr1	4232	4262	.	.	+
",
    "chr1	4233	4263	.	.	+
",
    "chr1	4234	4264	.	.	+
",
    "chr1	4235	4265	.	.	+
",
    "chr1	4237	4267	.	.	+
",
    "chr1	4237	4267	.	.	+
",
    "chr1	4239	4269	.	.	+
",
    "chr1	4239	4269	.	.	+
",
    "chr1	4244	4274	.	.	+
",
    "chr1	4245	4275	.	.	+
",
    "chr1	4247	4277	.	.	+
",
    "chr1	4247	4277	.	.	+
",
    "chr1	4248	4278	.	.	+
",
    "chr1	4249	4279	.	.	+
",
    "chr1	4249	4279	.	.	+
",
    "chr1	4250	4280	.	.	+
",
    "chr1	4251	4281	.	.	+
",
    "chr1	4252	4282	.	.	+
",
    "chr1	4253	4283	.	.	+
",
    "chr1	4254	4284	.	.	+
",
    "chr1	4255	4285	.	.	+
",
    "chr1	4257	4287	.	.	+
",
    "chr1	4258	4288	.	.	+
",
    "chr1	4259	4289	.	.	+
",
    "chr1	4259	4289	.	.	+
",
    "chr1	4260	4290	.	.	+
",
    "chr1	4260	4290	.	.	+
",
    "chr1	4260	4290	.	.	+
",
);

/// 60 single-end control reads, deliberately denser between the clusters.
const CTRL: &str = concat!(
    "chr1	760	790	.	.	+
",
    "chr1	897	927	.	.	+
",
    "chr1	1140	1170	.	.	+
",
    "chr1	1163	1193	.	.	+
",
    "chr1	1719	1749	.	.	+
",
    "chr1	1823	1853	.	.	+
",
    "chr1	1953	1983	.	.	+
",
    "chr1	2444	2474	.	.	+
",
    "chr1	2549	2579	.	.	+
",
    "chr1	3008	3038	.	.	+
",
    "chr1	3033	3063	.	.	+
",
    "chr1	3059	3089	.	.	+
",
    "chr1	3095	3125	.	.	+
",
    "chr1	3208	3238	.	.	+
",
    "chr1	3214	3244	.	.	+
",
    "chr1	3227	3257	.	.	+
",
    "chr1	3233	3263	.	.	+
",
    "chr1	3233	3263	.	.	+
",
    "chr1	3237	3267	.	.	+
",
    "chr1	3238	3268	.	.	+
",
    "chr1	3243	3273	.	.	+
",
    "chr1	3253	3283	.	.	+
",
    "chr1	3285	3315	.	.	+
",
    "chr1	3295	3325	.	.	+
",
    "chr1	3305	3335	.	.	+
",
    "chr1	3307	3337	.	.	+
",
    "chr1	3311	3341	.	.	+
",
    "chr1	3316	3346	.	.	+
",
    "chr1	3326	3356	.	.	+
",
    "chr1	3344	3374	.	.	+
",
    "chr1	3361	3391	.	.	+
",
    "chr1	3363	3393	.	.	+
",
    "chr1	3366	3396	.	.	+
",
    "chr1	3387	3417	.	.	+
",
    "chr1	3825	3855	.	.	+
",
    "chr1	4232	4262	.	.	+
",
    "chr1	4463	4493	.	.	+
",
    "chr1	4870	4900	.	.	+
",
    "chr1	5030	5060	.	.	+
",
    "chr1	5072	5102	.	.	+
",
    "chr1	5086	5116	.	.	+
",
    "chr1	5150	5180	.	.	+
",
    "chr1	5215	5245	.	.	+
",
    "chr1	5232	5262	.	.	+
",
    "chr1	5232	5262	.	.	+
",
    "chr1	5249	5279	.	.	+
",
    "chr1	5267	5297	.	.	+
",
    "chr1	5317	5347	.	.	+
",
    "chr1	5365	5395	.	.	+
",
    "chr1	5420	5450	.	.	+
",
    "chr1	5432	5462	.	.	+
",
    "chr1	5436	5466	.	.	+
",
    "chr1	5472	5502	.	.	+
",
    "chr1	5491	5521	.	.	+
",
    "chr1	5499	5529	.	.	+
",
    "chr1	5634	5664	.	.	+
",
    "chr1	6119	6149	.	.	+
",
    "chr1	6145	6175	.	.	+
",
    "chr1	6399	6429	.	.	+
",
    "chr1	6407	6437	.	.	+
",
);

/// `*_summits.bed` for `-p 0.01` (with control).
const EXPECTED_SUMMITS_P: &str = concat!(
    "chr1\t2044\t2045\tx_peak_1\t17.2525\n",
    "chr1\t2730\t2731\tx_peak_2\t8.64252\n",
    "chr1\t4316\t4317\tx_peak_3\t33.0656\n",
);

/// `*_peaks.narrowPeak` for `-p 0.01` (with control).
const EXPECTED_NARROW_P: &str = concat!(
    "chr1\t2008\t2179\tx_peak_1\t172\t.\t10.0782\t17.2525\t15.7352\t36\n",
    "chr1\t2611\t2769\tx_peak_2\t86\t.\t5.57868\t8.64252\t7.59137\t119\n",
    "chr1\t4209\t4403\tx_peak_3\t330\t.\t15.5315\t33.0656\t29.4211\t107\n",
);

/// `*_summits.bed` for `--broad -p 0.01` (with control).
const EXPECTED_BROAD_P: &str = concat!(
    "chr1\t2004\t2184\tx_peak_1\t137\t.\t8.40676\t13.7942\t12.4497\n",
    "chr1\t2604\t2772\tx_peak_2\t80\t.\t5.49343\t8.09344\t7.01856\n",
    "chr1\t4205\t4407\tx_peak_3\t221\t.\t11.23\t22.1983\t20.25\n",
);

/// `*_summits.bed` for `--call-summits -p 0.01` -- the smoothed sub-peaks, which sit
/// at different positions from the plain summit and carry the p-value score.
const EXPECTED_SUBPEAK_P: &str = concat!(
    "chr1\t2093\t2094\tx_peak_1\t17.2525\n",
    "chr1\t2691\t2692\tx_peak_2\t10.7513\n",
    "chr1\t4305\t4306\tx_peak_3\t31.3661\n",
);

/// `*_summits.bed` for `-q 0.01`: **two** peaks, scoring with the q-value column.
const EXPECTED_SUMMITS_Q: &str = concat!(
    "chr1\t2044\t2045\tx_peak_1\t15.7352\n",
    "chr1\t4316\t4317\tx_peak_2\t29.4211\n",
);

fn tmpdir(tag: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("macs3rs-pcut-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn write_fixture(dir: &std::path::Path) -> (String, String) {
    let treat = dir.join("treat.bed");
    let ctrl = dir.join("ctrl.bed");
    std::fs::write(&treat, TREAT).unwrap();
    std::fs::write(&ctrl, CTRL).unwrap();
    (
        treat.to_str().unwrap().to_string(),
        ctrl.to_str().unwrap().to_string(),
    )
}

/// Run `callpeak ... -n x` into a fresh directory and return the body of `file`
/// with the `#`-prefixed header lines removed (the reference and this port
/// necessarily disagree on the paths echoed in `# Command line:`).
fn callpeak(tag: &str, extra: &[&str]) -> (String, String, String) {
    let dir = tmpdir(tag);
    let (treat, ctrl) = write_fixture(&dir);
    let out = dir.join("out");
    let mut args: Vec<String> = vec![
        "callpeak".into(),
        "-t".into(),
        treat,
        "-c".into(),
        ctrl,
        "-g".into(),
        "20000".into(),
        "--nomodel".into(),
        "--extsize".into(),
        "150".into(),
        "-n".into(),
        "x".into(),
    ];
    args.extend(extra.iter().map(|s| s.to_string()));
    args.push("--outdir".into());
    args.push(out.to_str().unwrap().to_string());
    let st = std::process::Command::new(RUST)
        .args(&args)
        .output()
        .unwrap();
    assert!(
        st.status.success(),
        "callpeak {:?} failed: {}",
        extra,
        String::from_utf8_lossy(&st.stderr)
    );
    let read = |name: &str| match std::fs::read_to_string(out.join(name)) {
        Ok(text) => text
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| format!("{l}\n"))
            .collect::<String>(),
        // a narrow run writes no broadPeak, and a broad run no narrowPeak
        Err(_) => String::new(),
    };
    (
        read("x_summits.bed"),
        read("x_peaks.narrowPeak"),
        read("x_peaks.broadPeak"),
    )
}

#[test]
fn narrow_pvalue_cutoff_is_a_pvalue_cutoff() {
    let (summits, narrow, _) = callpeak("narrow-p", &["-p", "0.01"]);
    assert_eq!(summits, EXPECTED_SUMMITS_P);
    assert_eq!(narrow, EXPECTED_NARROW_P);
}

#[test]
fn the_pvalue_path_is_not_the_qvalue_path() {
    // the control for this whole test: under `-q 0.01` the oracle calls *two* peaks,
    // so the three-peak `-p 0.01` answer above can only come from the p-score track.
    let (summits, _, _) = callpeak("narrow-q", &["-q", "0.01"]);
    assert_eq!(summits, EXPECTED_SUMMITS_Q);
    assert_ne!(summits, EXPECTED_SUMMITS_P);
}

#[test]
fn broad_pvalue_cutoff_is_a_pvalue_cutoff() {
    let (_, _, broad) = callpeak("broad-p", &["-p", "0.01", "--broad"]);
    assert_eq!(broad, EXPECTED_BROAD_P);
}

#[test]
fn call_summits_with_a_pvalue_cutoff() {
    let (summits, _, _) = callpeak("subpeak-p", &["-p", "0.01", "--call-summits"]);
    assert_eq!(summits, EXPECTED_SUBPEAK_P);
}
