#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""보고서 그림 재생성기 (2026-09-24 · 경로 A).

정본 파일엔 재생성기를 붙일 것 — 인라인으로 그리고 그림만 저장하면 감사할 수 없다.
입력은 전부 판정기가 쓴 `_RESULT_*.json` 이고, 이 스크립트는 **읽기만** 한다.
숫자를 새로 계산하지 않는다: 판정 결과에 없는 값을 그리면 그림과 본문의 정본이 갈린다.

  F2  침입 Δ 대 μ — draw-first · find-first · draw-once (10-03 이름 바꿈) · 둘째 시조   (_RESULT_inv9 · _RESULT_inv10)
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

축 라벨 규약 (2026-09-29 · 투고 점검에서 정함 — 새 그림도 이걸 따른다)
  · 오류율 축은 한 문구로: "nominal error rate, $\\mu$ (%)"  (10-03: μ 는 글자당이 아니라 추첨당 값)
  · Δ · ΔΔ 계열은 무단위가 아니다 — "(ln units)" · "(ln units per tick)" 를 붙인다
  · 고리 길이 단위는 "ticks per copy attempt" 하나로(글자당 아님 — 헛손질이 있으면 갈린다)
  · 변수는 mathtext 로 이탤릭($\\mu$ · $\\Delta$ · $I$ · $R$ · $F$ · $U$ · $a$), 코드 이름은 $\\mathtt{...}$
  · 로그 눈금은 라벨이나 캡션이 말한다 · 천 단위는 쉼표(1,200)
  · 주사위 은유는 표 4 에만 — 축은 "error draws per letter written, $R$"
  · 선 종류(실선/점선) 같은 범례는 축 라벨에 넣지 않는다(범례·캡션)
  · **패널 안에는 통계도 설명도 두지 않는다** — 저널 표본 11장에서 패널 안 통계 0건(본문·캡션이 든다).
      남겨도 되는 것은 점의 이름(코드 이름)뿐이다.
  · 통계 수치를 그림 안에 또 적지 않는다 — 캡션과 자릿수가 어긋난다(09-29 에 F7c · FS1 에서 뺐다)
  · **범례·라벨이 데이터를 덮으면 안 된다** — 자리를 손으로 정하지 말고 `legend_clear(ax, ...)` 를 쓸 것
      (후보 자리를 다 그려 보고 겹치는 점을 세어 고르고, 어디에도 자리가 없으면 위쪽에 여유를 준다.
       09-29 에 손으로 정한 자리 여섯이 데이터를 덮고 있었다 — 사용자가 둘을 먼저 잡아냈다.)
  · 두 패널이 같은 계열을 쓰면 **그림 전체 범례 하나**로 둔다(09-29 사용자 결정) · 자리는 **그림 아래**(`loc="outside lower center"` · 10-03 사용자 — 패널 문자보다 위에 두지 않는다).
      패널 밖이라 데이터를 가릴 수 없고, (b) 만 봐도 색을 안다. F4 가 그 예다.
  · 한 그림 안에서 **같은 범례를 두 번 싣지 않는다**(F4(b) 가 (a) 와 같은 셋을 또 실어 y 축 글자에 붙었다 — 09-29 사용자 지적).
      캡션이 "colours as in (a)" 로 잇는다.
  · 범례가 커서 어디에도 안 들어가면 자리를 옮기지 말고 **항목을 줄인다**(캡션이 이미 말하는 항목은 범례에서 뺀다 — FS2).
  · 같은 양을 그린 패널끼리는 y 범위를 맞출 것(F8 은 범례 여유로 (a) 가 늘어나 (b) 와 어긋났다)
  · 라벨이 그림 밖으로 나가면 안 된다: 아래 한 줄로 기계 점검
      python3 - <<'EOF' (fm.save 를 가로채 fig.findobj(Text) 의 window_extent 를 fig.bbox 와 대조)

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
from matplotlib.lines import Line2D   # noqa: E402
from matplotlib.patches import Patch   # noqa: E402

D = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(os.path.dirname(D), "petri")   # 09-28 폴더 정리: 결과·스크립트·동결 해시는 petri/ (lab101 ~/petri/ 미러)
OUT = os.path.join(D, "figures")
os.makedirs(OUT, exist_ok=True)
PROV = []
MUS = [0.0, 0.001, 0.003, 0.005, 0.01]
COL = {"rascld": "#c0392b", "rsacld": "#2980b9", "racldx": "#7f8c8d", "racld": "#27ae60"}
LBL = {"rascld": r"$\mathtt{rascld}$ (loop 2)", "rsacld": r"$\mathtt{rsacld}$ (loop 3)",
       "racldx": r"$\mathtt{racldx}$ (loop 4, +no-op)", "racld": r"$\mathtt{racld}$ (loop 4)"}

plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 200, "savefig.bbox": "tight", "axes.linewidth": 0.7,
                     "xtick.major.width": 0.7, "ytick.major.width": 0.7,
                     "figure.constrained_layout.use": True,
                     "figure.constrained_layout.w_pad": 0.08})


LOCS = ("upper right", "upper left", "lower right", "lower left", "center right", "center left", "upper center", "lower center")


def _data_px(ax):
    """축 안의 그려진 것들을 화면 좌표 점으로. 기준선(점선)은 뺀다.

    09-29: 처음엔 선의 점만 셌다 — 막대 그림(4b · 7b · S2b)이 검사를 그냥 통과했다(사용자 지적).
    막대는 네 꼭짓점과 윗변 가운데를, 산점 묶음은 각 점을 넣는다.
    """
    pts = []
    for ln in ax.lines:
        if ln.get_linestyle() == ":" or len(ln.get_xdata()) == 0:
            continue
        pts += [tuple(q) for q in ax.transData.transform(list(zip(ln.get_xdata(), ln.get_ydata())))]
    for pa in ax.patches:                      # 막대
        bb = pa.get_window_extent()
        pts += [(bb.x0, bb.y0), (bb.x1, bb.y0), (bb.x0, bb.y1), (bb.x1, bb.y1),
                ((bb.x0 + bb.x1) / 2, bb.y1), ((bb.x0 + bb.x1) / 2, bb.y0)]
    for co in ax.collections:                  # 오차막대 · 산점 묶음
        try:
            off = co.get_offsets()
        except Exception:
            continue
        if off is not None and len(off):
            pts += [tuple(q) for q in ax.transData.transform(off)]
    return pts


def _mk(ax, handles, labels, loc, kw):
    if handles is None:
        return ax.legend(loc=loc, **kw)
    if labels is None:
        return ax.legend(handles=handles, loc=loc, **kw)
    return ax.legend(handles, labels, loc=loc, **kw)


def legend_clear(ax, handles=None, labels=None, locs=LOCS, **kw):
    """데이터를 가장 적게 가리는 자리에 범례를 놓는다.

    09-29: 자리를 손으로 정하다 여섯 군데에서 데이터를 덮었다(사용자가 둘을 잡아냈다).
    후보 자리를 다 그려 보고 겹치는 점 수를 세어 고른다 — 자료가 바뀌어도 따라온다.
    """
    fig = ax.figure

    def try_all():
        fig.canvas.draw()
        pts = _data_px(ax)
        best, best_hit = None, None
        for loc in locs:
            lg = _mk(ax, handles, labels, loc, kw)
            fig.canvas.draw()
            bb = lg.get_window_extent()
            hit = sum(1 for x, y in pts if bb.x0 - 1 <= x <= bb.x1 + 1 and bb.y0 - 1 <= y <= bb.y1 + 1)
            if best_hit is None or hit < best_hit:
                best, best_hit = loc, hit
            if hit == 0:
                break
        return best, best_hit

    best, hit = try_all()
    for _ in range(2):            # 어디에도 자리가 없으면 위쪽에 여유를 준다
        if not hit:
            break
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi * (hi / lo) ** 0.18 if ax.get_yscale() == "log" else hi + (hi - lo) * 0.18)
        best, hit = try_all()
    return _mk(ax, handles, labels, best, kw)


def load(tag):
    p = os.path.join(BASE, "_RESULT_%s.json" % tag)
    if not os.path.exists(p):
        sys.exit("없음: " + p)
    return json.load(open(p))


def note(fig, src, path, val):
    PROV.append("%-4s %-16s %-44s %s" % (fig, src, path, val))


def panel_labels(fig):
    """10-02(사용자 결정): "(a) 제목" 꼴 가운데 제목을 굵은 패널 문자 a + 왼쪽 정렬 짧은 제목으로 바꾼다.
    10-03(사용자 지적 셋 — 제목이 축 라벨보다 작다 · 제목이 너무 왼쪽 · 문자가 y 축 라벨과 겹친다):
      · 크기: 문자 9 pt 굵게 · 제목 8 pt(= 축 라벨, 지침 §2)
      · 문자는 축선에서 시작하고 제목이 4 pt 뒤에 붙는다 — 한 덩어리(사용자 "같이" · "아직 왼쪽")
      · 높이는 그림 안 모든 패널에 같은 값 — 제목 줄이 맞는다. 문자가 y 축 라벨과 겹치면 그림 전체를 함께 올린다
    그림마다 고치지 않고 여기 한 곳에서."""
    items = []
    for ax in fig.axes:
        m = re.match(r"^\((\w)\)\s*(.*)$", ax.get_title())
        if not m:
            continue
        ax.set_title("")
        items.append((ax, m.group(1), m.group(2)))
    if not items:
        return
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    pt = fig.dpi / 72.0                                          # 1 pt 의 화면 픽셀
    arts = []
    for ax, letter, title in items:
        lab = ax.annotate(letter, xy=(0, 1), xycoords="axes fraction", xytext=(0, 6), textcoords="offset points",
                          fontsize=9, fontweight="bold", ha="left", va="baseline", annotation_clip=False)
        tt = ax.annotate(title, xy=(0, 1), xycoords="axes fraction", xytext=(0, 6), textcoords="offset points",
                         fontsize=8, ha="left", va="baseline", annotation_clip=False) if title else None
        arts.append((ax, lab, tt))

    def place(dy):
        fig.canvas.draw()
        for ax, lab, tt in arts:
            ab = ax.get_window_extent(r)
            yb = ax.yaxis.get_tightbbox(r)                       # y 축이 없는 패널은 None
            lw = lab.get_window_extent(r).width / pt
            left_pt = ((yb.x0 - ab.x0) / pt) if yb is not None else 0.0
            x = 0.0                                               # 10-03(사용자: 아직 왼쪽): 문자가 축선에서 시작하고 제목이 4 pt 뒤에 붙는다
            lab.xyann = (x, dy)
            if tt is not None:
                tt.xyann = (lw + 4, dy)

    dy = 6.0
    place(dy)
    for _ in range(3):                                           # 문자–y 축 라벨 겹침이 없어질 때까지 그림 전체를 함께 올린다
        fig.canvas.draw()
        up = 0.0
        for ax, lab, tt in arts:
            yb = ax.yaxis.get_tightbbox(r)
            lb = lab.get_window_extent(r)
            if yb is not None and lb.y0 < yb.y1 and lb.x0 < yb.x1 and lb.x1 > yb.x0:
                up = max(up, (yb.y1 - lb.y0) / pt + 2)
        if up <= 0:
            break
        dy += up
        place(dy)


def save(fig, name):
    panel_labels(fig)
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    fig.savefig(os.path.join(OUT, name + ".png"))
    fig.savefig(os.path.join(OUT, name + ".svg"))   # 09-29: PDF 조판에 벡터로 넣으려고 추가
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
    rules = [("orig", z10, "draw-first"), ("ff", z10, "find-first"), ("mem", z9, "draw-once"), ("anc", z10, r"draw-first, $\mathtt{racld}$-founded backgrounds")]
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
        ax.plot(x, ys, "-", marker="os^D"[i], ms=3, lw=1.2,   # 10-02: 흑백에서도 갈리게 마커를 규칙마다 다르게(사용자 결정)
                color=["#c0392b", "#2980b9", "#8e44ad", "#e67e22"][i], label=lab)
    ax.axhline(0, color="k", lw=0.6, ls=":")
    ax.set_xlabel(r"nominal error rate, $\mu$ (%)")
    ax.set_ylabel(r"invasion $\Delta$ of $\mathtt{rascld}$ (ln units)")   # 정의는 캡션에 있다
    legend_clear(ax, frameon=False, fontsize=6.5)
    save(fig, "F2_invasion_by_rule")


# ---------------------------------------------------------------- F3
def f3():
    z = load("spec8")
    fig, axes = plt.subplots(1, 2, figsize=(6.0, 2.5), sharey=True)
    for ax, rule, title in ((axes[0], "orig", "(a) draw-first"), (axes[1], "ff", "(b) find-first")):
        for code in ("rascld", "rsacld", "racld"):   # 10-02: 범례를 고리 길이 순(2 · 3 · 4)으로
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
        ax.set_yticks([1, 3, 10]); ax.set_yticklabels(["1", "3", "10"])   # 10-02: 지수 표기 대신 보통 숫자
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        legend_clear(ax, frameon=False, fontsize=6.5) if rule == "ff" else None
        ax.set_xlabel("material density (letters per cell)")
    axes[0].set_ylabel(r"inflation, $I$")   # 정의 · 로그 눈금은 캡션이 말한다
    save(fig, "F3_inflation_by_density")


# ---------------------------------------------------------------- F4
def f4():
    z = load("an11")
    cells = z["A_cells"]
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.5))
    ax = axes[0]
    for code in ("rascld", "rsacld", "racld"):
        for d, ls in ((8, "-"), (14, "--")):
            xs, ys = [], []
            for a in (300, 600, 1200):
                c = cells["%s|%d|%d" % (code, d, a)]
                xs.append(a); ys.append(c["s"])
                note("F4", "an11", "A_cells/%s|%d|%d/s" % (code, d, a), "%.2f" % c["s"])
            ax.plot(xs, ys, ls, marker="o", ms=3, lw=1.2, color=COL[code],
                    label=LBL[code] if d == 8 else None)
    ax.set_yscale("log"); ax.set_yticks([1, 3, 10, 30]); ax.set_yticklabels(["1", "3", "10", "30"])   # 10-02: 지수 표기 대신 보통 숫자 · 코드 순서는 고리 길이 순
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xticks([300, 600, 1200]); ax.set_xticklabels(["300", "600", "1,200"])
    ax.set_xlabel(r"lifespan constant $a$ (ticks)")
    ax.set_ylabel(r"stalls per copied letter, $h$")
    ax.set_title("(a) stalls and lifespan", fontsize=7)
    # 09-29: 두 패널이 같은 세 코드를 쓰므로 범례는 그림 전체에 하나 둔다(아래 fig.legend · 사용자 결정). 10-03: 위치는 그림 아래로(사용자).
    ax = axes[1]
    names = [r"change in ln $h$" "\n(total)", "stall\ntime share", "org-ticks\nper birth", "letters\nper birth"]
    for j, code in enumerate(("rascld", "rsacld", "racld")):
        c = z["A_contrast|%s|8" % code]
        vals = [c[0][0], c[1][0], c[2][0], -c[3][0]]
        for i, v in enumerate(vals):
            ax.bar(i + (j - 1) * 0.26, v, width=0.25, color=COL[code], label=LBL[code] if i == 0 else None)
        note("F4", "an11", "A_contrast|%s|8" % code, "Δln s %.3f · share %.3f · tpb %.3f · lpb %.3f" % (c[0][0], c[1][0], c[2][0], c[3][0]))
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(range(4)); ax.set_xticklabels(names)   # 10-03: (a) 와 같은 크기
    ax.set_ylabel(r"change in ln $h$")
    ax.set_title("(b) exact decomposition", fontsize=7)
    # 09-29: (a) 와 같은 세 코드라 범례를 두 번 싣지 않는다 — 색은 캡션이 "(a) 와 같다" 로 잇는다.
    ax.set_ylim(top=max(z["A_contrast|%s|8" % q][2][0] for q in ("rascld", "rsacld", "racld")) * 1.5)
    fig.legend(handles=[Line2D([], [], marker="o", ls="-", ms=3, lw=1.2, color=COL[c], label=LBL[c])
                        for c in ("rascld", "rsacld", "racld")],
               loc="outside lower center", ncol=3, frameon=False, fontsize=7,   # 10-03(사용자): 패널 문자 위가 아니라 그림 아래
               handlelength=1.4, columnspacing=1.6, borderaxespad=0.2)
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
        ax.annotate(r"$\mathtt{%s}$" % c, (L, v[0]), textcoords="offset points", xytext=(5, 4), fontsize=6)
        note("F5", "an11", "C_delta|" + c, "%.3f [%.3f, %.3f]" % (v[0], v[1][0], v[1][1]))
    for c, mk in (("racldn", "s"), ("racldr", "D")):
        v = z["C_delta|" + c]
        ax.errorbar([4], [v[0]], yerr=[[v[0] - v[1][0]], [v[1][1] - v[0]]], fmt=mk, ms=4,
                    mfc="none", color="#8e44ad", capsize=2, lw=1)
        ax.annotate(r"$\mathtt{%s}$" % c, (4, v[0]), textcoords="offset points", xytext=(5, -8), fontsize=6)
        note("F5", "an11", "C_delta|" + c, "%.3f [%.3f, %.3f]" % (v[0], v[1][0], v[1][1]))
    ax.axhline(0, color="k", lw=0.6, ls=":")
    ax.set_xticks([2, 3, 4, 5]); ax.set_xlim(1.6, 5.6); ax.set_xlabel("copy-loop length (ticks per attempt)")
    ax.set_ylabel(r"invasion $\Delta$ at $\mu$ = 0 (ln units)")
    ax.set_title("(a) loop-length series", fontsize=7)
    ax.set_ylim(-1.75, 0.75)
    v = z["C_delta|racldr"]
    # 09-29: racldr 의 값(−3.52)과 축 밖이라는 사실은 캡션에 있다 — 화살표만 남긴다.
    ax.annotate("", (4.0, -1.70), xytext=(4.0, -0.95), arrowprops=dict(arrowstyle="->", color="#8e44ad", lw=0.8))
    ax.annotate(r"$\mathtt{racldr}$", (4.08, -1.32), ha="left", va="center", fontsize=6, color="#8e44ad")
    # 09-29: 화살표는 racldr 이 축 아래(−3.52)로 벗어났다는 표시다. 가로 위치는 고리 길이 4 여야 한다 —
    #        4.55 에 두면 가로축이 고리 길이인 그림에서 '고리 4.5' 로 읽힌다(사용자 지적).
    ax = axes[1]
    tab = b_table()
    # 10-02(사용자 결정): 직선 적합 점선을 뺐다 — 굽은 자료 위에 떠서 '모형이 못 맞춘다' 로 읽혔다. 기울기 · 비는 본문 3.7 이 말한다.
    #        색은 순서 있는 양이라 한 색조(회색)의 연한 → 진한으로 — 빨강 · 주황은 Fig. 2 에서 규칙의 색이다.
    for a, col in ((300, "#b3b3b3"), (600, "#737373"), (1200, "#252525")):
        xs = [m * 100 for m in MUS]
        ys = [tab[a][m] for m in MUS]
        for m in MUS:
            note("F5", "an11.txt", "B table lifespan %d · mu %g" % (a, m), "%+.3f" % tab[a][m])
        ax.plot(xs, ys, "o-", ms=3, lw=1.2, color=col, label="lifespan %s" % format(a, ","))
    ax.axhline(0, color="k", lw=0.6, ls=":")
    ax.set_xlabel(r"nominal error rate, $\mu$ (%)")
    ax.set_ylabel(r"invasion $\Delta$ (ln units)", fontsize=7)
    legend_clear(ax, frameon=False, fontsize=6.5)
    ax.set_title("(b) three lifespans", fontsize=7)
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
    ax.axvspan(0.05, 0.35, color="#f1c40f", alpha=0.10, lw=0)   # 띠의 뜻은 캡션(보충 S7)이 말한다
    ax.axhline(0, color="k", lw=0.6, ls=":")
    ax.set_xlabel(r"nominal error rate, $\mu$ (%)")
    ax.set_ylabel(r"invasion $\Delta$ (ln units)", fontsize=7)
    legend_clear(ax, frameon=False, fontsize=5.8, labelspacing=0.35)   # 실선·파선의 뜻도 캡션에 있다 — 코드 셋만 싣는다
    ax.set_title("(a) measurement and model", fontsize=7)
    ax = axes[1]
    codes = ["rascld", "rsacld", "racld", "racldx"]   # 10-02: 고리 길이 순
    for i, c in enumerate(codes):
        f = z["flow|" + c]
        ax.bar(i, f["U"], color=COL[c], width=0.6)
        note("FS2", "model12", "flow|%s/U" % c, "%.4f" % f["U"])
    ax.set_xticks(range(len(codes))); ax.set_xticklabels([r"$\mathtt{%s}$" % c for c in codes], fontsize=6.5)
    ax.set_ylabel(r"variant offspring per birth", fontsize=7)
    ax.set_title("(b) variant flow", fontsize=7)
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
        yt.append((y, r"$\mathtt{%s}$ (loop %d→%d)" % (b, L, L + 1)))
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
    ax.set_xlabel(r"$\Delta\Delta$ per added loop tick (ln units)")
    ax.set_xlim(-1.45, 0.25)
    legend_clear(ax, handles=[
        Line2D([], [], marker="o", ls="", ms=2.6, color="#34495e", label="reread"),
        Line2D([], [], marker="o", ls="", ms=2.6, color="#c0392b", label="pre-registered"),
        Patch(facecolor="#bdc3c7", alpha=0.25, label="pre-registered range"),
        Line2D([], [], ls="--", lw=0.8, color="#27ae60", label="series slope (Fig. 3a)"),   # 10-02: 패널 안 수치 뺌
    ], frameon=False, fontsize=5.2, handlelength=1.5, borderaxespad=0.25, labelspacing=0.35)
    ax.set_title("(a) matched pairs", fontsize=7)
    ax = axes[1]
    ax.axvline(0.15, color="#c0392b", lw=0.8, ls="--")
    # 10-02(사용자 결정): 잡음 바닥 띠 — 같은 코드 두 팔이 다를 수 있는 폭의 95% 분위(해부 23 · 첫 배경 세트 S1).
    #        `acld` 는 해부 23 C0 로 판정 밖이라 띠가 없다(캡션이 말한다).
    z23 = load("noise23")
    for i, (b, L, src) in enumerate(ORDER):
        y = len(ORDER) - 1 - i
        r = z23["rows"].get("S1|%s" % b)
        if r and r["judged"]:
            ax.barh(y, r["q95"], height=0.84, left=0, color="#bdc3c7", alpha=0.45, lw=0, zorder=0)
            note("F6", "noise23", "rows/S1|%s/q95" % b, "%.3f [%.3f, %.3f]" % (r["q95"], r["q95_ci"][0], r["q95_ci"][1]))
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
    ax.set_xlabel(r"$|\Delta\Delta|$, equal-loop pairs (ln units)")
    ax.set_xlim(-0.02, 0.73)
    legend_clear(ax, handles=[
        Line2D([], [], marker="o", ls="", ms=2.6, color="#34495e", label=r"inserted $\mathtt{n}$"),
        Line2D([], [], marker="o", ls="", ms=2.6, mfc="none", color="#34495e", label=r"inserted $\mathtt{x}$"),
        Line2D([], [], ls="--", lw=0.8, color="#c0392b", label="tolerance 0.15"),
        Patch(facecolor="#bdc3c7", alpha=0.45, label="noise floor"),
    ], frameon=False, fontsize=5.2, handlelength=1.5, borderaxespad=0.25, labelspacing=0.35)
    ax.set_title("(b) position control", fontsize=7)
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
    ax.set_xticks([2, 3, 4, 5]); ax.set_xlim(1.6, 5.9)
    ax.set_xlabel("copy-loop length (ticks per attempt)")
    ax.set_ylabel(r"error draws per letter written, $R$")
    legend_clear(ax, handles=[Line2D([], [], marker="o", ls="-", ms=2.8, lw=0.8, color=BC[b], label=r"$\mathtt{%s}$" % b)
                              for b, _, _ in ORDER], frameon=False, fontsize=5.2, ncol=2,
                 handlelength=1.2, columnspacing=0.8, labelspacing=0.3, borderaxespad=0.25)
    a1 = [v for k, v in Z14.items() if k.startswith("A1|")]
    a2 = Z14.get("A2_abs", {})
    # 09-29: 범위 둘(+0.16~+0.36 · ≤ 0.050)은 캡션에 있다 — 패널 안에서 또 말하지 않는다.
    note("F6", "pairs14", "A1|* (%d등급 쌍)" % len(a1), "ln R 차 %+.3f ~ %+.3f" % (min(v["d"] for v in a1), max(v["d"] for v in a1)))
    note("F6", "pairs14", "A2_abs.max", "%.4f" % a2.get("max", float("nan")))
    ax.set_title("(c) draws per letter", fontsize=7)
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
        ax.plot([m * 100 for m in MU], mean, "o-", color=BC[b], ms=2.6, lw=1.1, label=r"$\mathtt{%s}$" % b)   # 10-02: 짝 수 · 실패 수는 캡션으로(패널 안 수치 금지)
        note("F7", "pairs15", "B1|%s|* (%d짝)" % (b, len(rows)), "ΔΔ(μ) 평균 %s · 통과 %d/%d" % (" ".join("%+.2f" % x for x in mean), sum(1 for v in rows if v["pass"]), len(rows)))
    ax.set_xlabel(r"nominal error rate, $\mu$ (%)")
    ax.set_ylabel(r"$\Delta\Delta$, longer loop minus shorter (ln units)")
    ax.set_xticks([0, 0.1, 0.3, 0.5, 1.0]); ax.set_ylim(-1.05, 0.35)
    legend_clear(ax, fontsize=5.2, frameon=False, ncol=2, handlelength=1.4)
    ax.set_title("(a) error-rate series", fontsize=7)
    ax = axes[1]
    bars = [("draw-first", Z["C2|orig"], "#c0392b"), ("find-first", Z["C3|ff"], "#2980b9"), ("draw-once", Z["C2|rem"], "#8e44ad")]   # 10-02: Fig. 2 의 규칙 색 · 이름
    for i, (lab, v, col) in enumerate(bars):
        ax.bar(i, v["m"], color=col, width=0.62, alpha=0.9)
        ax.errorbar([i], [v["m"]], yerr=[[v["m"] - v["ci"][0]], [v["ci"][1] - v["m"]]], fmt="none", ecolor="k", capsize=2, lw=0.8)
        note("F7", "pairs15", "C2|orig · C3|ff · C2|rem" if i == 0 else "", "%s Δμ %+.3f [%+.3f, %+.3f]" % (lab, v["m"], v["ci"][0], v["ci"][1]))
    c1 = [("find-first", Z["C1|ff"]), ("draw-once", Z["C1|rem"])]
    for i, (lab, v) in enumerate(c1, start=1):
        ax.plot([i], [v["m"]], "D", ms=3.2, color="white", mec="k", mew=0.7, zorder=4)
        note("F7", "pairs15", "C1|%s" % ("ff" if i == 1 else "rem"), "μ0 차 %+.3f [%+.3f, %+.3f]" % (v["m"], v["ci"][0], v["ci"][1]))
    ax.axhline(0, color="k", lw=0.6, ls=":")
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels([b[0] for b in bars], fontsize=6.4)
    ax.set_ylabel(r"rise of $\Delta\Delta$, $\mu$ = 0 to 0.3% (ln units)")
    legend_clear(ax, handles=[
        Patch(facecolor="#bdc3c7", label=r"rise of $\Delta\Delta$ with error"),
        Line2D([], [], marker="D", ls="", ms=4, mfc="white", mec="k", mew=0.7,
               label=r"change of the $\mu$ = 0 difference"),
    ], frameon=False, fontsize=5.2, handlelength=1.2, borderaxespad=0.25, labelspacing=0.35)
    ax.set_title("(b) three copying rules", fontsize=7)
    # (c) 해부 18 — 순수 racld 공동체 · 밀도 8 대 32: 충실도 읽기(이웃 등급 ln R 차)와 μ 상승, 속도 읽기(|ΔΔ(0)|)
    ax = axes[2]
    d1 = sorted((k.split("|")[1], v) for k, v in Z18.items() if k.startswith("D1|"))
    x8 = [v["d8"] for _, v in d1]; x32 = [v["d32"] for _, v in d1]
    for i, (b, v) in enumerate(d1):
        ax.plot([0, 1], [v["d8"], v["d32"]], "-o", color=BC.get(b, "#34495e"), ms=2.8, lw=0.9, alpha=0.9)
        ax.annotate(r"$\mathtt{%s}$" % b, (1.03, v["d32"]), fontsize=5.0, color=BC.get(b, "#34495e"), va="center")
        note("F7", "pairs18", "D1|%s (d8 · d32 · ratio)" % b, "%+.3f · %+.3f · %.2f" % (v["d8"], v["d32"], v["ratio"] if v["ratio"] is not None else float("nan")))
    ax.axhline(0, color="k", lw=0.6, ls=":")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["density 8", "density 32"], fontsize=6.4); ax.set_xlim(-0.3, 1.55)
    ax.set_ylabel("ln ratio of draws, shorter / longer")
    D2 = Z18["D2"]; D3 = Z18["D3"]; D1b = Z18["D1b"]
    # 09-29: 통계 넷은 그림 안이 아니라 캡션으로 옮겼다(본문과 같은 자릿수로).
    note("F7", "pairs18", "D1b.median · D2.m8 · D2.m32 · D3.abs8 · D3.abs32", "%.3f · %+.3f · %+.3f · %.3f · %.3f" % (D1b["median"], D2["m8"], D2["m32"], D3["abs8"], D3["abs32"]))
    ax.set_title("(c) two densities", fontsize=7)
    save(fig, "F7_ladder_and_rules")


def f9():
    """진화 궤적 — 해부 17 (_RESULT_evo17). (a) 개체 가중 평균 고리 길이 L̄(t) · μ 네 층 · 접시 중앙값과 사분위 (b) 끝점 L̄ 대 μ · 접시마다."""
    Z = load("evo17")
    MU = [0.001, 0.003, 0.005, 0.01]
    MC = {0.001: "#6ece58", 0.003: "#1f958b", 0.005: "#375a8c", 0.01: "#471063"}   # 10-02: μ 는 순서 있는 양 — viridis 4단계(연한 쪽 0.1%) · 규칙 색(빨강 · 보라)과 겹치지 않게
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.6), gridspec_kw={"width_ratios": [1.5, 1]})
    ax = axes[0]
    for mu in MU:
        c = Z["curve|%s" % mu]
        t = [x / 1000 for x in c["t"]]
        ax.fill_between(t, c["q1"], c["q3"], color=MC[mu], alpha=0.15, lw=0)
        ax.plot(t, c["median"], "-", color=MC[mu], lw=1.2, label=r"$\mu$ = %g%%" % (mu * 100))
        note("F8", "evo17", "curve|%s (median · q1 · q3)" % mu, "끝 %.3f · 최대 %.3f" % (c["median"][-1], max(c["median"])))
    ax.axhline(2, color="#7f8c8d", lw=0.5, ls=":"); ax.axhline(4, color="#7f8c8d", lw=0.5, ls=":")
    ax.set_xlabel("time (thousands of ticks)")
    ax.set_ylabel("mean copy-loop length\n(ticks per attempt)")
    ax.set_ylim(1.8, 4.5)
    legend_clear(ax, fontsize=5.6, frameon=False, ncol=2)
    ax.set_title("(a) loop length over time", fontsize=7)
    ax = axes[1]
    for i, mu in enumerate(MU):
        v = Z["V1|%s" % mu]["values"]
        xs = [i + (j - (len(v) - 1) / 2.0) * (0.5 / max(len(v), 1)) for j in range(len(v))]
        ax.plot(xs, v, "o", ms=2.2, color=MC[mu], alpha=0.75, lw=0)
        ax.plot([i - 0.32, i + 0.32], [Z["V1|%s" % mu]["median"]] * 2, "-", color="k", lw=1.0)
        note("F8", "evo17", "V1|%s" % mu, "끝점 L̄ 중앙 %.3f · 접시 %d" % (Z["V1|%s" % mu]["median"], Z["V1|%s" % mu]["n"]))
    ax.set_xticks(range(len(MU))); ax.set_xticklabels(["%g" % (mu * 100) for mu in MU]); ax.set_xlabel(r"nominal error rate, $\mu$ (%)")
    ax.set_ylabel("mean copy-loop length\n(ticks per attempt)"); ax.set_ylim(1.8, 4.5)
    # 09-29: Spearman 은 캡션으로 — 패널 안에 통계를 적지 않는다(저널 표본 11장에서 0건).
    note("F8", "evo17", "V1.rho · ci", "%.3f [%.3f, %.3f]" % (Z["V1"]["rho"], Z["V1"]["ci"][0], Z["V1"]["ci"][1]))
    ax.set_title("(b) endpoint by error rate", fontsize=7)
    axes[1].set_ylim(axes[0].get_ylim())   # 같은 양이므로 두 패널의 y 범위를 맞춘다(범례 여유로 (a) 가 늘어날 수 있다)
    save(fig, "F8_evolutionary_trajectory")


def fS1():
    """보충 그림 S1 — 해부 19 (_RESULT_pairs19 · 시험은 PETRI_FIG19=pilot19). 수명 사다리 × 밀도: (a) 합친 |ΔΔ(0)| (b) 성장 끝 빈 칸 몫."""
    tag = os.environ.get("PETRI_FIG19", "pairs19")
    if not os.path.exists(os.path.join(BASE, "_RESULT_%s.json" % tag)):
        print("  FS1 건너뜀 — _RESULT_%s.json 없음(판정 전)" % tag); return
    Z = load(tag)
    AGES = [150, 300, 600]
    DC = {"8": "#a6a6a6", "32": "#303030"}   # 10-02: 밀도는 순서 있는 양 — 회색 두 단계(초록 · 빨강은 코드 색이다)
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.5), gridspec_kw={"width_ratios": [1.25, 1]})
    ax = axes[0]
    for d in ("8", "32"):
        m = [Z["A|d%s|a%d" % (d, a)]["mean"] for a in AGES]
        lo = [Z["A|d%s|a%d" % (d, a)]["ci"][0] for a in AGES]; hi = [Z["A|d%s|a%d" % (d, a)]["ci"][1] for a in AGES]
        ax.errorbar(AGES, m, yerr=[[x - l for x, l in zip(m, lo)], [h - x for x, h in zip(m, hi)]], fmt="o-", color=DC[d], ms=3, lw=1.1, capsize=2, label="density %s" % d)
        note("FS1", tag, "A|d%s|a{150,300,600}.mean" % d, " ".join("%.3f" % x for x in m))
    ax.set_xscale("log"); ax.set_xticks(AGES); ax.set_xticklabels(["150", "300", "600"]); ax.minorticks_off()
    ax.set_xlabel(r"lifespan constant $a$ (ticks)")
    ax.set_ylabel(r"$|\Delta\Delta|$ at $\mu$ = 0 (ln units)")
    ax.set_ylim(bottom=0)
    v1, v2 = Z["V1"], Z["V2"]
    # 09-29: V1 · V2 수치는 캡션(보충 S6)이 세 자리로 들고 있어 그림에서 뺐다.
    note("FS1", tag, "V1.diff · V2.contrast", "%+.3f [%+.3f, %+.3f] · %+.3f [%+.3f, %+.3f]" % (v1["diff"], v1["ci"][0], v1["ci"][1], v2["contrast"], v2["ci"][0], v2["ci"][1]))
    legend_clear(ax, fontsize=5.8, frameon=False)
    ax.set_title("(a) lifespan series", fontsize=7)
    ax = axes[1]
    for d in ("8", "32"):
        v = [Z["cell|d%s|a%d" % (d, a)]["vac_median"] * 100 for a in AGES]
        ax.plot(AGES, v, "o-", color=DC[d], ms=3, lw=1.1, label="density %s" % d)
        note("FS1", tag, "cell|d%s|a{150,300,600}.vac_median" % d, " ".join("%.3f" % (x / 100) for x in v))
    ax.set_xscale("log"); ax.set_xticks(AGES); ax.set_xticklabels(["150", "300", "600"]); ax.minorticks_off()
    ax.set_xlabel(r"lifespan constant $a$ (ticks)")
    ax.set_ylabel("empty cells at the end of growth (%)"); ax.set_ylim(0, 60)
    p1 = Z["P1"]
    # 09-29: 빈 칸 비(1.79 · 0.57)도 캡션에 있다.
    note("FS1", tag, "P1.r150 · P1.r600", "%.3f · %.3f" % (p1["r150"], p1["r600"]))
    ax.set_title("(b) vacancies", fontsize=7)
    save(fig, "FS1_lifespan_ladder")


def fS3():
    """보충 그림 S3 — 해부 29 (_RESULT_evo29). 끝점 평균 고리 길이 대 명목 μ · 규칙 셋 · 접시마다 점 · 중앙값 가로선."""
    Z = load("evo29")
    RULES = [("DF", "draw-first", "#c0392b", [0.001, 0.003, 0.005, 0.01]),
             ("FF", "find-first", "#2980b9", [0.001, 0.003, 0.005, 0.01]),
             ("DO", "draw-once", "#8e44ad", [0.001, 0.003, 0.005, 0.01, 0.03, 0.05])]   # 규칙 색은 그림 2 와 같다
    MU = [0.001, 0.003, 0.005, 0.01, 0.03, 0.05]
    off = {"DF": -0.22, "FF": 0.0, "DO": 0.22}
    fig, ax = plt.subplots(figsize=(4.6, 2.7))
    for key, lab, col, mus in RULES:
        for mu in mus:
            c = Z["cell|%s|%s" % (key, mu)]
            i = MU.index(mu) + off[key]
            v = c["values"]
            xs = [i + (j - (len(v) - 1) / 2.0) * (0.16 / max(len(v), 1)) for j in range(len(v))]
            ax.plot(xs, v, "o", ms=1.8, color=col, alpha=0.6, lw=0)
            ax.plot([i - 0.09, i + 0.09], [c["median"]] * 2, "-", color="k", lw=0.9)
            note("FS3", "evo29", "cell|%s|%s" % (key, mu), "중앙 %.3f · 평균 %.3f · 접시 %d" % (c["median"], c["mean"], c["n"]))
        ax.plot([], [], "o", ms=3, color=col, label=lab)
    ax.axhline(2, color="#7f8c8d", lw=0.5, ls=":"); ax.axhline(4, color="#7f8c8d", lw=0.5, ls=":")
    ax.set_xticks(range(len(MU))); ax.set_xticklabels(["%g" % (m * 100) for m in MU])
    ax.set_xlabel(r"nominal error rate, $\mu$ (%)")
    ax.set_ylabel("mean copy-loop length at 50,000 ticks\n(ticks per attempt)")
    ax.set_ylim(1.8, 4.5)
    legend_clear(ax, fontsize=6, frameon=False)
    save(fig, "FS3_evolution_by_rule")


def fS4():
    """보충 그림 S4 — 갱신 계산(사후 · _RESULT_theory30 · 원자료 _RESULT_inv28). (a) q 별 Δ(μ) (b) Δ(μ) − Δ(0) 를 초과 추첨 공급 x = (52/33)Σ(q)μ 에 대해."""
    Z = load("theory30")
    Q = [0.0, 0.25, 0.5, 0.75, 1.0]
    QC = {0.0: "#fde725", 0.25: "#5ec962", 0.5: "#21918c", 0.75: "#3b528b", 1.0: "#440154"}   # q 는 순서 있는 양 — viridis 5단계
    pts = Z["collapse"]["points"]
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.6))
    ax = axes[0]
    for q in Q:
        ps = sorted([p for p in pts if p["q"] == q], key=lambda p: p["mu"])
        d0 = ps[0]["d0"]
        mus = [0.0] + [p["mu"] * 100 for p in ps]
        m = [d0[0]] + [p["d"][0] for p in ps]
        lo = [d0[1]] + [p["d"][1] for p in ps]; hi = [d0[2]] + [p["d"][2] for p in ps]
        ax.fill_between(mus, lo, hi, color=QC[q], alpha=0.18, lw=0)
        ax.plot(mus, m, "o-", color=QC[q], ms=2.6, lw=1.1, label=r"$q$ = %g" % q, mec="#333333" if q == 0.0 else QC[q], mew=0.4)
        note("FS4", "theory30", "collapse.points q=%g (d)" % q, " ".join("%.3f" % x for x in m))
    ax.axhline(0, color="k", lw=0.5, ls=":")
    ax.set_xlabel(r"nominal error rate, $\mu$ (%)")
    ax.set_ylabel(r"invasion $\Delta$ of $\mathtt{rascld}$ (ln units)")
    legend_clear(ax, fontsize=5.8, frameon=False)
    ax.set_title("(a) invasion by redraw probability", fontsize=7)
    ax = axes[1]
    for q in Q:
        ps = sorted([p for p in pts if p["q"] == q], key=lambda p: p["mu"])
        ax.plot([p["x"] for p in ps], [p["dd"] for p in ps], "o-", color=QC[q], ms=2.6, lw=0.9, mec="#333333" if q == 0.0 else QC[q], mew=0.4)
        note("FS4", "theory30", "collapse.points q=%g (x · dd)" % q, " ".join("%.4f/%+.3f" % (p["x"], p["dd"]) for p in ps))
    ax.axhline(0, color="k", lw=0.5, ls=":")
    ax.set_xlabel("estimated excess realised letter changes\nper offspring, founder over resident")
    ax.set_ylabel(r"$\Delta(\mu) - \Delta(0)$ (ln units)")
    ax.set_title("(b) against estimated excess letter changes", fontsize=7)
    axes[1].set_ylim(axes[0].get_ylim()[0] - 0.4, 0.15)
    save(fig, "FS4_renewal_collapse")


FIGS = {"F2": f2, "F3": f3, "F4": f4, "F5": f5, "F6": f7, "F7": f8, "F8": f9, "FS1": fS1, "FS2": f6, "FS3": fS3, "FS4": fS4}


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
    lines = PROV
    pp = os.path.join(OUT, "_PROVENANCE.txt")
    if only and os.path.exists(pp):
        # 10-02: --only 가 출처 파일을 그 그림 줄만 남기고 덮었다(138줄 소실 · 조용한 실패). 다른 그림의 줄은 그대로 둔다.
        keep = [ln for ln in io.open(pp, encoding="utf-8").read().split("\n")[2:] if ln.strip() and ln.split()[0] != only]
        order = {k: i for i, k in enumerate(FIGS)}
        lines = sorted(keep + PROV, key=lambda ln: order.get(ln.split()[0], 99))   # 안정 정렬 — 그림 안 순서는 유지
    io.open(pp, "w", encoding="utf-8").write(
        "그림에 그린 값의 출처 — 판정기 기록의 JSON 경로. 새로 계산한 값은 없다.\n\n" + "\n".join(lines) + "\n")
    PROV[:] = lines
    print("  출처 %d줄 → figures/_PROVENANCE.txt" % len(PROV))


if __name__ == "__main__":
    main()
