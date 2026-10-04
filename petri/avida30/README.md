# Experiment 30 — the copying rules in Avida

`material-draw.patch` applies to Avida (https://github.com/devosoft/avida) at commit 47f13da and adds one instruction, `h-copy-mat`, with conserved per-cell material and the three copying rules (draw-first, find-first, draw-once). With `MATERIAL_MODE 0` it behaves as `h-copy`. `apply_patch.py` regenerates the patch edits; `cfg/` holds the configuration used.

The patch modifies Avida, which is distributed under the GNU Lesser General Public License v3 or later; the patch and the files in this folder derived from Avida (`material-draw.patch`, `apply_patch.py`, `cfg/avida.cfg`) are released under the same licence (LGPL-3.0-or-later), not under the MIT licence of the rest of this repository.

Scripts: `a30_launch.sh` (checks the frozen hashes, branch, patch and binary before launching), `a30_e0.py` (validation), `a30_e1.py` (single-code census), `a30_e2.py` (invasion assays), `a30_analyze.py` (decision criteria), `a30lib.py`. Results: `_RESULT_a30_e*.txt` and the per-run tables `_RESULT_a30_e1_rows.json.gz`, `_RESULT_a30_e2_arms.json.gz`; pilots `_PILOT_a30_*`. Supplementary S28 of the manuscript describes the design and results.
