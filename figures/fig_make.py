#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""보고서 그림 재생성기 (2026-09-24 · 경로 A).

정본 파일엔 재생성기를 붙일 것 — 인라인으로 그리고 그림만 저장하면 감사할 수 없다.
입력은 전부 판정기가 쓴 `_RESULT_*.json` 이고, 이 스크립트는 **읽기만** 한다.
숫자를 새로 계산하지 않는다: 판정 결과에 없는 값을 그리면 그림과 본문의 정본이 갈린다.

  F2  침입 Δ 대 μ — roll-first · find-first · remember-roll · 둘째 시조   (_RESULT_inv9 · _RESULT_inv10)
  F3  부풀림 I = R × F — 밀도 × 규칙                                      (_RESULT_spec8)
  F4  헛손질/성공 대 수명과 분해                                          (_RESULT_an11)
  F5  (a) 고리 길이 사다리 (b) 수명별 침입 곡선                           (_RESULT_an11)
  F7  (a)(b) 짝지은 삽입 ΔΔ · 음성 대조 (c) 글자당 주사위 횟수 직접 계수     (_RESULT_pairs13 · pairs13_racldfam · pairs14)
  F8  (a) 짝 안 μ 사다리 ΔΔ(μ) (b) 세 규칙의 충실도 몫                        (_RESULT_pairs15)
  F9  (a) 진화 궤적 L̄(t) · μ 네 층 (b) 끝점 L̄ 대 μ                            (_RESULT_evo17)
  F6  해석 모형 — 예측 대 실측 · 변이 흐름                                (_RESULT_model12)
  FS1 보충 — 수명 사다리 × 밀도: |ΔΔ(0)| 와 빈 칸 몫                          (_RESULT_pairs19)
  F1 은 손 그림이라 여기서 만들지 않는다 — 소스는 `figures/F1_world_and_rules.svg`(09-26),
  PDF · PNG 는 `rsvg-convert -f pdf|png` 로 만든다. 판정 기록에서 나온 값이 아니므로 _PROVENANCE 에도 없다.

검산: 그린 값마다 출처 (파일 · JSON 경로)를 `figures/_PROVENANCE.txt` 에 적는다.
실행: python3 report/fig_make.py [--only F4]
"""
import io
import json
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

D = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(os.path.dirname(D), "petri")   # 09-28 폴더 정리: 결과·스크립트·동결 해시는 petri/ (lab101 ~/petri/ 미러)
OUT = os.path.join(D, "figures")
os.makedirs(OUT, exist_ok=True)
PROV = []
MUS = [0.0, 0.001, 0.003, 0.005, 0.01]
COL = {"rascld": "#c0392b", "rsacld": "#2980b9", "racldx": "#7f8c8d", "racld": "#27ae60"}
LBL = {"rascld": "rascld (loop 2)", "rsacld": "rsacld (loop 3)", "racldx": "racldx (loop 4, +no-op)", "racld": "racld (loop 4)"}

plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 200, "savefig.bbox": "tight", "axes.linewidth": 0.7,
                     "xtick.major.width": 0.7, "ytick.major.width": 0.7,
                     "figure.constrained_layout.use": True,
                     "figure.constrained_layout.w_pad": 0.08})


def load(tag):
    p = os.path.join(BASE, "_RESULT_%s.json" % tag)
    if not os.path.exists(p):
        sys.exit("없음: " + p)
    return json.load(open(p))


def note(fig, src, path, val):
    PROV.append("%-4s %-16s %-44s %s" % (fig, src, path, val))


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    fig.savefig(os.path.join(OUT, name + ".png"))
    plt.close(fig)
    print("  " + name)


def b_table():
    """해부 11 B 의 μ별 평균 — 판정기가 인쇄한 표에서 읽는다(JSON 에는 기울기·절편만 있다).
    여기서 새로 계산하지 않는다: 그림과 본문이 같은 값을 쓰게 하려면 판정기 출력이 정본이다."""
    txt = io.open(os.path.join(BASE, "_RESULT_an11.txt"), encoding="utf-8").read()
    head = [l for l in txt.split("\n") if l.strip().startswith("μ |")]
    if not head:
        sys.exit("_RESULT_an11.txt 에서 B 표를 못 찾았다")
    ages = [int(x) for x in re.findall(r"(\d+)", head[0])]
    out = {a: {} for a in ages}
    start = txt.index(head[0])
    for line in txt[start:].split("\n")[1:]:
        m = re.match(r"\s*([0-9.]+)\s*\|(.*)$", line)
        if not m:
            break
        mu = float(m.group(1))
        vals = [float(x) for x in re.findall(r"[+-][0-9.]+", m.group(2))]
        if len(vals) != len(ages):
            sys.exit("B 표의 칸 수가 수명 수와 다르다: %r" % line)
        for a, v in zip(ages, vals):
            out[a][mu] = v
    missing = [(a, mu) for a in ages for mu in MUS if mu not in out[a]]
    if missing:
        sys.exit("B 표에 없는 칸: %s" % missing[:4])
    return out


# ---------------------------------------------------------------- F2
def f2():
    z9, z10 = load("inv9"), load("inv10")
    rules = [("orig", z10, "roll-first"), ("ff", z10, "find-first"), ("mem", z9, "remember-roll"), ("anc", z10, "roll-first, racld-founded")]
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    for i, (r, z, lab) in enumerate(rules):
        ys, lo, hi = [], [], []
        for mu in MUS:
            k = "%s|rascld|%g" % (r, mu)
            v = z["delta"][k]
            ys.append(v[0]); lo.append(v[1]); hi.append(v[2])
            note("F2", "inv9/inv10", "delta/" + k, "%.3f [%.3f, %.3f]" % tuple(v[:3]))
        x = [m * 100 for m in MUS]
        ax.fill_between(x, lo, hi, color=["#c0392b", "#2980b9", "#8e44ad", "#e67e22"][i], alpha=0.12, lw=0)
        ax.plot(x, ys, "o-", ms=3, lw=1.2, color=["#c0392b", "#2980b9", "#8e44ad", "#e67e22"][i], label=lab)
    ax.axhline(0, color="k", lw=0.6, ls=":")
    ax.set_xlabel("nominal error rate μ per letter (%)")
    ax.set_ylabel("invasion Δ of rascld\n(log-ratio change over 3,000 ticks)")
    ax.legend(frameon=False, fontsize=6.5, loc="lower left")
    save(fig, "F2_invasion_by_rule")


# ---------------------------------------------------------------- F3
def f3():
    z = load("spec8")
    fig, axes = plt.subplots(1, 2, figsize=(6.0, 2.5), sharey=True)
    for ax, rule, title in ((axes[0], "orig", "roll-first"), (axes[1], "ff", "find-first")):
        for code in ("racld", "rascld", "rsacld"):
            xs, ys = [], []
            for d in (4, 8, 32):
                k = "%s/%s/%d/0.001" % (rule, code, d)
                if k not in z.get("cells", {}):
                    continue
                c = z["cells"][k]
                xs.append(d); ys.append(c["I"])
                note("F3", "spec8", "cells/" + k + "/I", "%.3f" % c["I"])
            if xs:
                ax.plot(xs, ys, "o-", ms=3, lw=1.2, color=COL[code], label=LBL[code])
        ax.axhline(1, color="k", lw=0.6, ls=":")
        ax.set_xscale("log"); ax.set_xticks([4, 8, 32]); ax.set_xticklabels(["4", "8", "32"])
        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax.set_yscale("log"); ax.set_title(title, fontsize=8)
        ax.legend(frameon=False, fontsize=6.5) if rule == "ff" else None
        ax.set_xlabel("material density (letters per cell)")
    axes[0].set_ylabel("inflation I = R × F\n(realised ÷ nominal variant supply)")
    save(fig, "F3_inflation_by_density")


# ---------------------------------------------------------------- F4
def f4():
    z = load("an11")
    cells = z["A_cells"]
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.5))
    ax = axes[0]
    for code in ("racld", "rascld", "rsacld"):
        for d, ls in ((8, "-"), (14, "--")):
            xs, ys = [], []
            for a in (300, 600, 1200):
                c = cells["%s|%d|%d" % (code, d, a)]
                xs.append(a); ys.append(c["s"])
                note("F4", "an11", "A_cells/%s|%d|%d/s" % (code, d, a), "%.2f" % c["s"])
            ax.plot(xs, ys, ls, marker="o", ms=3, lw=1.2, color=COL[code],
                    label=LBL[code] if d == 8 else None)
    ax.set_yscale("log"); ax.set_xticks([300, 600, 1200])
    ax.set_xlabel("lifespan constant (ticks)")
    ax.set_ylabel("stalls per copied letter")
    ax.legend(frameon=False, fontsize=6.2, loc="upper left", title="solid: density 8   dashed: 14", title_fontsize=6)
    ax = axes[1]
    names = ["Δln s\n(total)", "stall\ntime share", "org-ticks\nper birth", "letters\nper birth"]
    for j, code in enumerate(("racld", "rascld", "rsacld")):
        c = z["A_contrast|%s|8" % code]
        vals = [c[0][0], c[1][0], c[2][0], -c[3][0]]
        for i, v in enumerate(vals):
            ax.bar(i + (j - 1) * 0.26, v, width=0.25, color=COL[code], label=LBL[code] if i == 0 else None)
        note("F4", "an11", "A_contrast|%s|8" % code, "Δln s %.3f · share %.3f · tpb %.3f · lpb %.3f" % (c[0][0], c[1][0], c[2][0], c[3][0]))
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(range(4)); ax.set_xticklabels(names, fontsize=6)
    ax.set_ylabel("Δln, lifespan 1200 − 300 (density 8)", fontsize=7)
    ax.legend(frameon=False, fontsize=6.2, loc="upper center", ncol=3, columnspacing=0.7, handlelength=1.0)
    ax.set_ylim(top=max(z["A_contrast|%s|8" % q][2][0] for q in ("racld", "rascld", "rsacld")) * 1.5)
    save(fig, "F4_lifespan_budget")


# ---------------------------------------------------------------- F5
def f5():
    z = load("an11")
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.5))
    ax = axes[0]
    loop = {"rascld": 2, "rsacld": 3, "racldx": 4, "nracld": 5}
    for c, L in sorted(loop.items(), key=lambda kv: kv[1]):
        v = z["C_delta|" + c]
        ax.errorbar([L], [v[0]], yerr=[[v[0] - v[1][0]], [v[1][1] - v[0]]], fmt="o", ms=4,
                    color=COL.get(c, "#34495e"), capsize=2, lw=1)
        ax.annotate(c, (L, v[0]), textcoords="offset points", xytext=(5, 4), fontsize=6)
        note("F5", "an11", "C_delta|" + c, "%.3f [%.3f, %.3f]" % (v[0], v[1][0], v[1][1]))
    for c, mk in (("racldn", "s"), ("racldr", "D")):
        v = z["C_delta|" + c]
        ax.errorbar([4], [v[0]], yerr=[[v[0] - v[1][0]], [v[1][1] - v[0]]], fmt=mk, ms=4,
                    mfc="none", color="#8e44ad", capsize=2, lw=1)
        ax.annotate(c, (4, v[0]), textcoords="offset points", xytext=(5, -8), fontsize=6)
        note("F5", "an11", "C_delta|" + c, "%.3f [%.3f, %.3f]" % (v[0], v[1][0], v[1][1]))
    ax.axhline(0, color="k", lw=0.6, ls=":")
    ax.set_xticks([2, 3, 4, 5]); ax.set_xlim(1.6, 5.6); ax.set_xlabel("copy-loop length (ticks)")
    ax.set_ylabel("invasion Δ at μ = 0")
    ax.set_ylim(-1.75, 0.75)
    v = z["C_delta|racldr"]
    ax.annotate("racldr %.2f\n(binds a letter the\ncommunity already uses)" % v[0], (3.55, -1.62),
                ha="center", va="center", fontsize=5.8, color="#8e44ad")
    ax.annotate("", (4.55, -1.72), xytext=(4.55, -1.15), arrowprops=dict(arrowstyle="->", color="#8e44ad", lw=0.8))
    ax = axes[1]
    tab = b_table()
    for a, col in ((300, "#95a5a6"), (600, "#e67e22"), (1200, "#c0392b")):
        b = z["B"]["%d|rascld" % a]
        note("F5", "an11", "B/%d|rascld" % a, "mu0 %.3f · slope %.1f" % (b["mu0"][0], b["slope"][0]))
        xs = [m * 100 for m in MUS]
        ax.plot(xs, [b["mu0"][0] + b["slope"][0] * m for m in MUS], "--", lw=0.9, color=col, alpha=0.75)
        ys = [tab[a][m] for m in MUS]
        for m in MUS:
            note("F5", "an11.txt", "B table lifespan %d · mu %g" % (a, m), "%+.3f" % tab[a][m])
        ax.plot(xs, ys, "o-", ms=3, lw=1.2, color=col, label="lifespan %d" % a)
    ax.axhline(0, color="k", lw=0.6, ls=":")
    ax.set_xlabel("nominal error rate μ (%)")
    ax.set_ylabel("invasion Δ of rascld, racld-founded worlds\n(points: measured · dashed: fitted line)", fontsize=7)
    ax.legend(frameon=False, fontsize=6.5)
    save(fig, "F5_loop_and_lifespan")


# ---------------------------------------------------------------- F6
def f6():
    z = load("model12")
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.5))
    ax = axes[0]
    for code in ("rascld", "rsacld", "racldx"):
        me, pr, xs = [], [], []
        for mu in MUS[1:]:
            c = z["cell|%s|%g" % (code, mu)]
            xs.append(mu * 100); me.append(c["meas"]); pr.append(c["M1"])
            note("FS2", "model12", "cell|%s|%g" % (code, mu), "meas %.3f · M1 %.3f · hit %s" % (c["meas"], c["M1"], c["judged"] and c["hit"]))
        ax.plot(xs, me, "o-", ms=3, lw=1.2, color=COL[code], label=LBL[code])
        ax.plot(xs, pr, "^--", ms=3, lw=1.0, color=COL[code], alpha=0.65)
    ax.axvspan(0.05, 0.35, color="#f1c40f", alpha=0.10, lw=0)
    ax.axhline(0, color="k", lw=0.6, ls=":")
    ax.set_xlabel("nominal error rate μ (%)")
    ax.set_ylabel("invasion Δ  (solid: measured,\ndashed: parameter-free model M1)", fontsize=7)
    ax.legend(frameon=False, fontsize=6.2, loc="lower left")
    ax.set_title("shaded: pre-registered judging range", fontsize=6)
    ax = axes[1]
    codes = ["racld", "rascld", "rsacld", "racldx"]
    for i, c in enumerate(codes):
        f = z["flow|" + c]
        ax.bar(i, f["U"], color=COL[c], width=0.6)
        note("FS2", "model12", "flow|%s/U" % c, "%.4f" % f["U"])
    ax.set_xticks(range(len(codes))); ax.set_xticklabels(codes, fontsize=6.5)
    ax.set_ylabel("variant offspring per birth, U", fontsize=7)
    ax.set_title("measured inside the invasion, μ = 0.3%", fontsize=6.5)
    save(fig, "FS2_model_vs_measured")


ORDER = [("rascld", 2, "re-read"), ("rsacld", 3, "re-read"), ("racld", 4, "re-read"),
         ("racldx", 4, "re-read"), ("rascled", 2, "registered"), ("acld", 3, "registered")]


def f7():
    """짝지은 삽입 — 여섯 밑 코드. (a) 고리 한 틱의 ΔΔ · (b) 음성 대조 · (c) 글자당 주사위 횟수 R 을 직접 센 것(해부 14).
    세 기록에서 읽는다: 사전등록한 둘째 가족은 pairs13, `racld` 가족은 같은 판정기를 사후로 돌린 pairs13_racldfam,
    R 은 pairs14(invbud 계수 · μ 0 · 같은 배경 20)."""
    Z = {"registered": load("pairs13"), "re-read": load("pairs13_racldfam")}
    Z14 = load("pairs14")
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 2.6),
                             gridspec_kw={"width_ratios": [1.25, 1, 1.05]})
    ax = axes[0]
    ax.axvspan(-1.0, -0.2, color="#bdc3c7", alpha=0.25, lw=0)
    ax.axvline(0, color="k", lw=0.6, ls=":")
    ax.axvline(-0.510, color="#27ae60", lw=0.8, ls="--")
    yt = []
    for i, (b, L, src) in enumerate(ORDER):
        y = len(ORDER) - 1 - i
        yt.append((y, "%s (loop %d→%d)" % (b, L, L + 1)))
        col = "#c0392b" if src == "registered" else "#34495e"
        pr = [(k, v) for k, v in Z[src].items() if k.startswith("pair|%s|n|" % b)]
        assert pr, "짝 없음 " + b
        for j, (k, v) in enumerate(sorted(pr)):
            assert abs(v["per_tick"] - v["dd"]) < 1e-9, "틱 차가 1 이 아니다 " + k
            off = (j - (len(pr) - 1) / 2.0) * (0.62 / max(len(pr), 1))
            ax.errorbar([v["dd"]], [y + off], xerr=[[v["dd"] - v["ci"][0]], [v["ci"][1] - v["dd"]]],
                        fmt="o", ms=2.6, color=col, capsize=1.2, lw=0.7, alpha=0.95)
        d = [v["dd"] for _, v in pr]
        note("F6", "pairs13" if src == "registered" else "pairs13_racldfam",
             "pair|%s|n|* (%d짝)" % (b, len(pr)), "틱당 %+.3f ~ %+.3f" % (min(d), max(d)))
    for y in range(len(ORDER) - 1):
        ax.axhline(y + 0.5, color="#ecf0f1", lw=0.6, zorder=0)
    ax.set_yticks([y for y, _ in yt]); ax.set_yticklabels([t for _, t in yt], fontsize=6.4)
    ax.set_ylim(-0.6, len(ORDER) - 0.4)
    ax.set_xlabel("ΔΔ per added loop tick (letter `n` held fixed)")
    ax.set_xlim(-1.45, 0.25)
    ax.annotate("ladder −0.510", (-0.510, -0.55), ha="center", va="bottom",
                fontsize=5.6, color="#27ae60")
    ax.annotate("pre-registered range", (-1.43, len(ORDER) - 0.45), ha="left", va="top",
                fontsize=5.6, color="#7f8c8d")
    ax.set_title("(a) 60 matched pairs, all negative", fontsize=7)
    ax = axes[1]
    ax.axvline(0.15, color="#c0392b", lw=0.8, ls="--")
    for i, (b, L, src) in enumerate(ORDER):
        y = len(ORDER) - 1 - i
        for ch, mfc, dy in (("n", None, 0.16), ("x", "none", -0.16)):
            g = sorted(abs(v["dd"]) for k, v in Z[src].items() if k.startswith("neg|%s|%s|" % (b, ch)))
            if not g:
                continue
            ax.plot(g, [y + dy] * len(g), "o", ms=2.6, mfc=mfc,
                    color="#c0392b" if src == "registered" else "#34495e", lw=0, alpha=0.9)
            note("F6", "pairs13" if src == "registered" else "pairs13_racldfam",
                 "neg|%s|%s|* (%d쌍)" % (b, ch, len(g)),
                 "최대 |ΔΔ| %.3f · 0.15 밖 %d" % (max(g), sum(1 for x in g if x > 0.15)))
    for y in range(len(ORDER) - 1):
        ax.axhline(y + 0.5, color="#ecf0f1", lw=0.6, zorder=0)
    ax.set_yticks([len(ORDER) - 1 - i for i in range(len(ORDER))]); ax.set_yticklabels([])
    ax.set_ylim(-0.6, len(ORDER) - 0.4)
    ax.set_xlabel("|ΔΔ| between equal-loop pairs")
    ax.set_xlim(-0.02, 0.73)
    ax.annotate("tolerance 0.15", (0.155, -0.5), ha="left", va="bottom", fontsize=5.6, color="#c0392b")
    ax.annotate("filled `n`  ·  open `x`", (0.71, len(ORDER) - 0.5), ha="right", va="top", fontsize=5.6,
                color="#7f8c8d")
    ax.set_title("(b) the control that position is neutral", fontsize=7)
    # (c) 해부 14 — 같은 40 코드에서 글자당 주사위 횟수를 직접 센 것. 같은 고리 길이의 배치는 같은 x 에 겹친다.
    ax = axes[2]
    BC = {"rascld": "#c0392b", "rsacld": "#2980b9", "racld": "#27ae60", "racldx": "#7f8c8d", "rascled": "#8e44ad", "acld": "#d35400"}
    for i, (b, L0, src) in enumerate(ORDER):
        rows = [(k.split("|")[2], v) for k, v in Z14.items() if k.startswith("code|%s|" % b)]
        assert rows, "해부 14 코드 없음 " + b
        jit = (i - (len(ORDER) - 1) / 2.0) * 0.09
        cls = {}
        for c, v in rows:
            cls.setdefault(v["L"], []).append(v["R"])
            ax.plot([v["L"] + jit], [v["R"]], "o", ms=2.8, color=BC[b], alpha=0.85, lw=0, zorder=3)
        xs = sorted(cls)
        ax.plot([x + jit for x in xs], [sum(cls[x]) / len(cls[x]) for x in xs], "-", color=BC[b], lw=0.8, alpha=0.7, zorder=2)
        note("F6", "pairs14", "code|%s|* (%d코드)" % (b, len(rows)),
             "R %.2f ~ %.2f · L %s" % (min(v["R"] for _, v in rows), max(v["R"] for _, v in rows), xs))
        ax.annotate(b, (xs[-1] + jit + 0.06, sum(cls[xs[-1]]) / len(cls[xs[-1]])), fontsize=5.4, color=BC[b], va="center")
    ax.set_xticks([2, 3, 4, 5]); ax.set_xlim(1.6, 5.9)
    ax.set_xlabel("copy-loop length (ticks per attempt)")
    ax.set_ylabel("die draws per letter written, R")
    a1 = [v for k, v in Z14.items() if k.startswith("A1|")]
    a2 = Z14.get("A2_abs", {})
    ax.annotate("adjacent classes: ln-ratio %+.2f to %+.2f\nequal loop, other position: |ln-ratio| ≤ %.3f" %
                (min(v["d"] for v in a1), max(v["d"] for v in a1), a2.get("max", float("nan"))),
                (0.98, 0.97), xycoords="axes fraction", ha="right", va="top", fontsize=5.4, color="#7f8c8d")
    note("F6", "pairs14", "A1|* (%d등급 쌍)" % len(a1), "ln R 차 %+.3f ~ %+.3f" % (min(v["d"] for v in a1), max(v["d"] for v in a1)))
    note("F6", "pairs14", "A2_abs.max", "%.4f" % a2.get("max", float("nan")))
    ax.set_title("(c) draws per letter, counted directly", fontsize=7)
    save(fig, "F6_matched_pairs")


def f8():
    """짝 안의 μ 사다리와 세 규칙 — 해부 15 (_RESULT_pairs15). (a) 밑 코드마다 짝 평균 ΔΔ(μ) (b) 규칙별 합친 Δμ 와 μ 0 차."""
    Z = load("pairs15")
    Z18 = load("pairs18")
    MU = [0.0, 0.001, 0.003, 0.005, 0.01]
    BC = {"rascld": "#c0392b", "rsacld": "#2980b9", "racld": "#27ae60", "racldx": "#7f8c8d", "rascled": "#8e44ad", "acld": "#d35400"}
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 2.6), gridspec_kw={"width_ratios": [1.35, 1, 1]})
    ax = axes[0]
    ax.axhline(0, color="k", lw=0.6, ls=":")
    for b, L0, src in ORDER:
        rows = [v for k, v in Z.items() if k.startswith("B1|%s|" % b)]
        assert rows, "B1 없음 " + b
        for v in rows:
            ax.plot([m * 100 for m in MU], [v["dd_by_mu"][str(m)] if str(m) in v["dd_by_mu"] else v["dd_by_mu"][m] for m in MU],
                    "-", color=BC[b], lw=0.5, alpha=0.25)
        mean = [sum((v["dd_by_mu"].get(str(m), v["dd_by_mu"].get(m)) for v in rows)) / len(rows) for m in MU]
        ax.plot([m * 100 for m in MU], mean, "o-", color=BC[b], ms=2.6, lw=1.1, label="%s (%d pairs%s)" % (b, len(rows), "" if all(v["pass"] for v in rows) else ", %d fail" % sum(1 for v in rows if not v["pass"])))
        note("F7", "pairs15", "B1|%s|* (%d짝)" % (b, len(rows)), "ΔΔ(μ) 평균 %s · 통과 %d/%d" % (" ".join("%+.2f" % x for x in mean), sum(1 for v in rows if v["pass"]), len(rows)))
    ax.set_xlabel("nominal error rate μ (% per letter)"); ax.set_ylabel("ΔΔ, longer loop minus shorter")
    ax.set_xticks([0, 0.1, 0.3, 0.5, 1.0]); ax.set_ylim(-1.05, 0.35)
    ax.legend(fontsize=5.2, frameon=False, loc="upper left", ncol=2, handlelength=1.4)
    ax.set_title("(a) one loop tick across the error ladder, 60 pairs", fontsize=7)
    ax = axes[1]
    bars = [("original", Z["C2|orig"], "#34495e"), ("find-first", Z["C3|ff"], "#7f8c8d"), ("remember-roll", Z["C2|rem"], "#c0392b")]
    for i, (lab, v, col) in enumerate(bars):
        ax.bar(i, v["m"], color=col, width=0.62, alpha=0.9)
        ax.errorbar([i], [v["m"]], yerr=[[v["m"] - v["ci"][0]], [v["ci"][1] - v["m"]]], fmt="none", ecolor="k", capsize=2, lw=0.8)
        note("F7", "pairs15", "C2|orig · C3|ff · C2|rem" if i == 0 else "", "%s Δμ %+.3f [%+.3f, %+.3f]" % (lab, v["m"], v["ci"][0], v["ci"][1]))
    c1 = [("find-first", Z["C1|ff"]), ("remember-roll", Z["C1|rem"])]
    for i, (lab, v) in enumerate(c1, start=1):
        ax.plot([i], [v["m"]], "D", ms=3.2, color="white", mec="k", mew=0.7, zorder=4)
        note("F7", "pairs15", "C1|%s" % ("ff" if i == 1 else "rem"), "μ0 차 %+.3f [%+.3f, %+.3f]" % (v["m"], v["ci"][0], v["ci"][1]))
    ax.axhline(0, color="k", lw=0.6, ls=":")
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels([b[0] for b in bars], fontsize=6.4)
    ax.set_ylabel("rise of ΔΔ from μ = 0 to 0.3%, pooled")
    ax.annotate("bars: fidelity reading\ndiamonds: change of the μ = 0 difference", (0.98, 0.97), xycoords="axes fraction", ha="right", va="top", fontsize=5.4, color="#7f8c8d")
    ax.set_title("(b) the rule decides which reading survives", fontsize=7)
    # (c) 해부 18 — 순수 racld 공동체 · 밀도 8 대 32: 충실도 읽기(이웃 등급 ln R 차)와 μ 상승, 속도 읽기(|ΔΔ(0)|)
    ax = axes[2]
    d1 = sorted((k.split("|")[1], v) for k, v in Z18.items() if k.startswith("D1|"))
    x8 = [v["d8"] for _, v in d1]; x32 = [v["d32"] for _, v in d1]
    for i, (b, v) in enumerate(d1):
        ax.plot([0, 1], [v["d8"], v["d32"]], "-o", color=BC.get(b, "#34495e"), ms=2.8, lw=0.9, alpha=0.9)
        ax.annotate(b, (1.03, v["d32"]), fontsize=5.0, color=BC.get(b, "#34495e"), va="center")
        note("F7", "pairs18", "D1|%s (d8 · d32 · ratio)" % b, "%+.3f · %+.3f · %.2f" % (v["d8"], v["d32"], v["ratio"] if v["ratio"] is not None else float("nan")))
    ax.axhline(0, color="k", lw=0.6, ls=":")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["density 8", "density 32"], fontsize=6.4); ax.set_xlim(-0.3, 1.55)
    ax.set_ylabel("ln-ratio of draws, shorter over longer loop")
    D2 = Z18["D2"]; D3 = Z18["D3"]; D1b = Z18["D1b"]
    ax.annotate("pure resident community\nmedian draws per letter %.2f → %.2f\nrise of ΔΔ to 0.3%%: %+.2f → %+.2f\n|ΔΔ| at μ = 0: %.2f → %.2f" % (
        D1b.get("median8", float("nan")), D1b["median"], D2["m8"], D2["m32"], D3["abs8"], D3["abs32"]),
        (0.03, 0.62), xycoords="axes fraction", ha="left", va="center", fontsize=5.2, color="#7f8c8d")
    note("F7", "pairs18", "D1b.median · D2.m8 · D2.m32 · D3.abs8 · D3.abs32", "%.3f · %+.3f · %+.3f · %.3f · %.3f" % (D1b["median"], D2["m8"], D2["m32"], D3["abs8"], D3["abs32"]))
    ax.set_title("(c) material made plentiful", fontsize=7)
    save(fig, "F7_ladder_and_rules")


def f9():
    """진화 궤적 — 해부 17 (_RESULT_evo17). (a) 개체 가중 평균 고리 길이 L̄(t) · μ 네 층 · 접시 중앙값과 사분위 (b) 끝점 L̄ 대 μ · 접시마다."""
    Z = load("evo17")
    MU = [0.001, 0.003, 0.005, 0.01]
    MC = {0.001: "#2980b9", 0.003: "#8e44ad", 0.005: "#d35400", 0.01: "#c0392b"}
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.6), gridspec_kw={"width_ratios": [1.5, 1]})
    ax = axes[0]
    for mu in MU:
        c = Z["curve|%s" % mu]
        t = [x / 1000 for x in c["t"]]
        ax.fill_between(t, c["q1"], c["q3"], color=MC[mu], alpha=0.15, lw=0)
        ax.plot(t, c["median"], "-", color=MC[mu], lw=1.2, label="μ = %g%%" % (mu * 100))
        note("F8", "evo17", "curve|%s (median · q1 · q3)" % mu, "끝 %.3f · 최대 %.3f" % (c["median"][-1], max(c["median"])))
    ax.axhline(2, color="#7f8c8d", lw=0.5, ls=":"); ax.axhline(4, color="#7f8c8d", lw=0.5, ls=":")
    ax.set_xlabel("ticks (thousands)"); ax.set_ylabel("mean copy-loop length of the population")
    ax.set_ylim(1.8, 4.5); ax.legend(fontsize=5.6, frameon=False, loc="lower right", ncol=2)
    ax.set_title("(a) from one two-tick founder, twenty dishes per rate", fontsize=7)
    ax = axes[1]
    for i, mu in enumerate(MU):
        v = Z["V1|%s" % mu]["values"]
        xs = [i + (j - (len(v) - 1) / 2.0) * (0.5 / max(len(v), 1)) for j in range(len(v))]
        ax.plot(xs, v, "o", ms=2.2, color=MC[mu], alpha=0.75, lw=0)
        ax.plot([i - 0.32, i + 0.32], [Z["V1|%s" % mu]["median"]] * 2, "-", color="k", lw=1.0)
        note("F8", "evo17", "V1|%s" % mu, "끝점 L̄ 중앙 %.3f · 접시 %d" % (Z["V1|%s" % mu]["median"], Z["V1|%s" % mu]["n"]))
    ax.set_xticks(range(len(MU))); ax.set_xticklabels(["%g" % (mu * 100) for mu in MU]); ax.set_xlabel("nominal error rate μ (% per letter)")
    ax.set_ylabel("mean loop length at 50,000 ticks"); ax.set_ylim(1.8, 4.5)
    ax.annotate("Spearman %.2f [%.2f, %.2f]" % (Z["V1"]["rho"], Z["V1"]["ci"][0], Z["V1"]["ci"][1]), (0.03, 0.97), xycoords="axes fraction", ha="left", va="top", fontsize=5.6, color="#7f8c8d")
    note("F8", "evo17", "V1.rho · ci", "%.3f [%.3f, %.3f]" % (Z["V1"]["rho"], Z["V1"]["ci"][0], Z["V1"]["ci"][1]))
    ax.set_title("(b) where the population ends up", fontsize=7)
    save(fig, "F8_evolutionary_trajectory")


def fS1():
    """보충 그림 S1 — 해부 19 (_RESULT_pairs19 · 시험은 PETRI_FIG19=pilot19). 수명 사다리 × 밀도: (a) 합친 |ΔΔ(0)| (b) 성장 끝 빈 칸 몫."""
    tag = os.environ.get("PETRI_FIG19", "pairs19")
    if not os.path.exists(os.path.join(BASE, "_RESULT_%s.json" % tag)):
        print("  FS1 건너뜀 — _RESULT_%s.json 없음(판정 전)" % tag); return
    Z = load(tag)
    AGES = [150, 300, 600]
    DC = {"8": "#27ae60", "32": "#c0392b"}
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.5), gridspec_kw={"width_ratios": [1.25, 1]})
    ax = axes[0]
    for d in ("8", "32"):
        m = [Z["A|d%s|a%d" % (d, a)]["mean"] for a in AGES]
        lo = [Z["A|d%s|a%d" % (d, a)]["ci"][0] for a in AGES]; hi = [Z["A|d%s|a%d" % (d, a)]["ci"][1] for a in AGES]
        ax.errorbar(AGES, m, yerr=[[x - l for x, l in zip(m, lo)], [h - x for x, h in zip(m, hi)]], fmt="o-", color=DC[d], ms=3, lw=1.1, capsize=2, label="density %s" % d)
        note("FS1", tag, "A|d%s|a{150,300,600}.mean" % d, " ".join("%.3f" % x for x in m))
    ax.set_xscale("log"); ax.set_xticks(AGES); ax.set_xticklabels(["150", "300", "600"]); ax.minorticks_off()
    ax.set_xlabel("lifespan a (ticks; lifetime = a + U[0, a))"); ax.set_ylabel("|ΔΔ| at μ = 0, pooled over 60 pairs")
    ax.set_ylim(bottom=0)
    v1, v2 = Z["V1"], Z["V2"]
    ax.annotate("density 32, 150 − 300: %+.2f [%+.2f, %+.2f]\ncontrast vs density 8: %+.2f [%+.2f, %+.2f]" % (v1["diff"], v1["ci"][0], v1["ci"][1], v2["contrast"], v2["ci"][0], v2["ci"][1]),
                (0.03, 0.04), xycoords="axes fraction", ha="left", va="bottom", fontsize=5.4, color="#7f8c8d")
    note("FS1", tag, "V1.diff · V2.contrast", "%+.3f [%+.3f, %+.3f] · %+.3f [%+.3f, %+.3f]" % (v1["diff"], v1["ci"][0], v1["ci"][1], v2["contrast"], v2["ci"][0], v2["ci"][1]))
    ax.legend(fontsize=5.8, frameon=False, loc="upper left")
    ax.set_title("(a) the speed reading along a lifespan ladder", fontsize=7)
    ax = axes[1]
    for d in ("8", "32"):
        v = [Z["cell|d%s|a%d" % (d, a)]["vac_median"] * 100 for a in AGES]
        ax.plot(AGES, v, "o-", color=DC[d], ms=3, lw=1.1, label="density %s" % d)
        note("FS1", tag, "cell|d%s|a{150,300,600}.vac_median" % d, " ".join("%.3f" % (x / 100) for x in v))
    ax.set_xscale("log"); ax.set_xticks(AGES); ax.set_xticklabels(["150", "300", "600"]); ax.minorticks_off()
    ax.set_xlabel("lifespan a (ticks)"); ax.set_ylabel("empty cells at the end of growth (%)"); ax.set_ylim(0, 60)
    p1 = Z["P1"]
    ax.annotate("density 32 vacancy ratio\n150/300: %.2f · 600/300: %.2f" % (p1["r150"], p1["r600"]), (0.97, 0.55), xycoords="axes fraction", ha="right", va="center", fontsize=5.4, color="#7f8c8d")
    note("FS1", tag, "P1.r150 · P1.r600", "%.3f · %.3f" % (p1["r150"], p1["r600"]))
    ax.set_title("(b) the manipulation: vacancies", fontsize=7)
    save(fig, "FS1_lifespan_ladder")


FIGS = {"F2": f2, "F3": f3, "F4": f4, "F5": f5, "F6": f7, "F7": f8, "F8": f9, "FS1": fS1, "FS2": f6}


def main():
    only = None
    if "--only" in sys.argv:
        only = sys.argv[sys.argv.index("--only") + 1]
    for k, fn in FIGS.items():
        if only and k != only:
            continue
        try:
            fn()
        except KeyError as e:
            print("  🔴 %s — 판정 기록에 없는 키 %s (판정기가 그 값을 안 쓴다 → 본문에도 쓰면 안 된다)" % (k, e))
    io.open(os.path.join(OUT, "_PROVENANCE.txt"), "w", encoding="utf-8").write(
        "그림에 그린 값의 출처 — 판정기 기록의 JSON 경로. 새로 계산한 값은 없다.\n\n" + "\n".join(PROV) + "\n")
    print("  출처 %d줄 → figures/_PROVENANCE.txt" % len(PROV))


if __name__ == "__main__":
    main()
