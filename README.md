# Code Petri Dish

Simulator, instruments, launch and judging scripts, frozen hashes and result files for the manuscript
*When a stalled copy draws again: copy-loop length couples replication speed to realised mutation rate under conserved material* (MinSeo Kim and Jae-Ho Shin, 2026).

## What is here

| Folder | Contents |
|---|---|
| `petri/sim.js` | The world (v0.3.6). One file, no dependencies. The browser trial and the server runs use the same file. v0.3.6 adds one option, `redrawP` (experiment 28); with its default the trajectory is identical to v0.3.5. |
| `petri/mini26.py`, `petri/mini27.py` | Two minimal models outside the simulator, in Python with NumPy (Supplementary S19): a well-mixed global pool, and local material on a 40 × 25 grid. `petri/mini31_*` (experiment 31) repeats the local model over a series reaching a nominal 4%, importing `mini27.py` and `mini27_analyze.py` unchanged. |
| `petri/p4win_*`, `petri/galA_window.py` | Post hoc recount of predictions P4a and P4b of the two minimal models (Tables S10 and S11; Supplementary S3) with draws compared at the same time in the assay. `galA_window.py` imports `mini27.py` unchanged and adds counters kept in 50-tick bins; the injected arms rerun identically to the original runs. Definitions were fixed and hashed (`p4win_frozen.sha256`) before the rerun; results in `_RESULT_p4win.*`, the weighting sensitivity in `_RESULT_p4win_posthoc.txt`. |
| `petri/avida30/` | Experiment 30: the three copying rules added to Avida as a patch (`material-draw.patch`, LGPL-3.0-or-later as a modification of Avida; see the folder's README), its validation, launch and analysis scripts, and results. |
| `petri/theory30_renewal.py` | The renewal account of the crossings (Supplementary S27), computed post hoc from existing result files into `_RESULT_theory30.*`. |
| `petri/*.js` | Instruments: `replay.js` (run one dish; `--loops 1` census of loop lengths), `dms.js` (inject one code at a time and follow the marked lineage), `invbud.js` (count copy outcomes per lineage: draws, stalls, deletions), `lineage.js`, `budget.js`, `found.js`, `mono.js`, `region.js`, `lifehist.js`. |
| `petri/*_launch.sh` | Launchers. Each writes a job list and runs the instrument over it with `xargs -P`. |
| `petri/*_analyze.py`, `petri/*_check.py` | Judging scripts and regression checks. Every number in the manuscript is written by one of these into `_RESULT_*.json` / `.txt`. |
| `petri/*_frozen.sha256` | SHA-256 of the simulator, launcher and judge frozen before each judging run; amended lines are appended, never replaced. |
| `petri/_RESULT_*` | Result files, one pair per judging run. |
| `figures/` | `fig_make.py` regenerates every plotted figure from the result files; `_PROVENANCE.txt` lists the JSON path behind every plotted value. The rendered PNG and PDF files are not deposited here from v1.4 on, because the journal of submission asks that the figure files of an article not be deposited in a public repository, and that material deposited here not be duplicated in the electronic supplementary material; run the script to obtain them. `F1_world_and_rules.svg` is kept because Figure 1 is hand-drawn and no script regenerates it. |

## Reproducing

Requires Node ≥ 20 and Python 3 (NumPy for the minimal models, matplotlib for figures).

```
cd petri
node sim.js --help                      # world options
node replay.js --seed 3 --mu 0.01       # one dish
node dms.js --cond mat --seed 3 --wt racld --arms list --codes nracld,rnacld --exact 1
python3 pairs19_analyze.py .            # a judge; reads _RESULT_* written by the launcher's jobs
python3 ../figures/fig_make.py          # all figures
```

The launchers assume a Linux host with many cores (the judging runs used 60 parallel jobs); paths inside them refer to `~/petri`.
Backgrounds are grown deterministically from the seed, so a run is reproducible byte for byte on the same Node version (v20.20.2 was used).

## Pre-registration

Each judging run was preceded by a note fixing the cells, the guards, the judging lines and the analysis script; the launcher and judge were hashed (`*_frozen.sha256`) before launch and re-checked after. Pilots ran on seeds disjoint from the judging seeds. The notes are in the project's research vault and are summarised in the manuscript's Methods 5.6 and Supplementary Table S2; frozen copies of the notes for experiments 23–27 are included (`prereg*_frozen.*`). For experiments 28–31 the hash list covers the vault note itself, which is not deposited and was later extended with the results. Outputs of the minimal models, like dish records, are not deposited; they regenerate from the seeds in the result files (NumPy 2.4.3 was used).

## Versions

- v1.0: experiments up to 25.
- v1.1: adds experiments 26–28 (the two minimal models and the re-draw probability series), `sim.js` v0.3.6, and the revised figures.
- v1.2: adds experiment 29 (evolution under draw-once and find-first; `replay.js --remember-die`), experiment 30 (the copying rules in Avida), the renewal account (Figure S3), the hash file of experiment 13, the revised title and figures.
- v1.3: after an independent code review — Figure S3 (evolution under the three rules; the theory figure becomes S4), independent recounts (`evo29_recount_independent.py`, `_RECOUNT29_independent.txt`, `avida30/_REVIEW30_*`), notes appended to the frozen hash files of experiments 29–30, and a corrected docstring in `theory30_renewal.py`. No result changed.
- v1.4: removes the rendered figure files (PNG and PDF), for the reason given in the `figures/` row above. `fig_make.py` and `_PROVENANCE.txt` are unchanged, so every plotted figure regenerates from the result files, and `F1_world_and_rules.svg` is retained as the drawing source of Figure 1. No result changed. The earlier tags here, and the Zenodo records of v1.0–v1.3, still contain the rendered files.
- v1.5: adds the post hoc recount of P4a and P4b in the two minimal models at the same time in the assay (`p4win_*`, `galA_window.py`, `_RESULT_p4win.*`). The pre-specified verdicts in `_RESULT_mini26.*` and `_RESULT_mini27.*` are unchanged; the manuscript reports both. Also adds experiment 31 (the local minimal model over a series reaching a nominal 4%; `mini31_*`, `dms31_frozen.sha256`, `_RESULT_mini31.*`), cited in Supplementary S19 but missing from v1.4. Supplementary section numbers in this README updated.

## License

MIT (see `LICENSE`), except `petri/avida30/material-draw.patch`, `apply_patch.py` and `cfg/avida.cfg`, which modify or derive from Avida and are LGPL-3.0-or-later.
