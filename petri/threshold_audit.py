# -*- coding: utf-8 -*-
"""문턱 민감도 감사 — 🟡 사후(판정 아님). 설계검토-2026-09-17.md J5 처방 · 사용자 문턱 검토용.
Wilson 범주 판정마다 '접시 하나(k ± 1)가 바뀌면 범주가 뒤집히나' 를 센다. 자료는 저장된 판정 출력(json)에서 읽는다 — 새 실행 없음.

  사건 문턱(Q1 · _RESULT_q1fix.json): 범주 = 매번(하한 ≥ 0.90) · 거의 안(상한 ≤ 0.10) · 우연 — 두 묶음(1~100 · 101~200) 각각
  되감기 2 Q2′ 정체(_RESULT_rewind2.json q2_id): 안정성 중앙 < 0.90 떠돈다 · 하한 ≥ 0.90 한 곳 · 상한 ≤ 0.50 갈라짐 · 그 밖 대체로 한 곳
  되감기 2 Q3′ 길(q3_path): 하한 ≥ 0.90 매번 · 상한 ≤ 0.10 드묾 · 그 밖 우연
  해부 8 A(_RESULT_rewind8.json): 재료 먼저 racld 몫 Wilson 상한 < 0.5 · 씨앗 공유 몫 Q2′ 규칙
실행: python3 threshold_audit.py  (볼트 폴더에서)  →  표준 출력 · _RESULT_threshold_audit.txt 로 저장한다
"""
import json
import math


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return (max(0.0, c - h), min(1.0, c + h))


def cat_q1(k, n):
    lo, hi = wilson(k, n)
    return "매번" if lo >= 0.90 else ("거의 안" if hi <= 0.10 else "우연")


def cat_q2(k, n):
    lo, hi = wilson(k, n)
    return "한 곳" if lo >= 0.90 else ("갈라짐" if hi <= 0.50 else "대체로 한 곳")


def cat_q3(k, n):
    lo, hi = wilson(k, n)
    return "매번" if lo >= 0.90 else ("드묾" if hi <= 0.10 else "우연")


rows = []


def audit(label, k, n, fn):
    c = fn(k, n)
    cm = fn(max(0, k - 1), n) if k > 0 else c
    cp = fn(min(n, k + 1), n) if k < n else c
    flip = (cm != c) or (cp != c)
    lo, hi = wilson(k, n)
    rows.append((label, k, n, lo, hi, c, cm, cp, flip))


q1 = json.load(open("_RESULT_q1fix.json"))
for cond in ("space", "energy", "mat"):
    for ev, d in q1["conds"][cond].items():
        if not all(x in d for x in ("k1", "n1", "k2", "n2")):
            continue
        audit("Q1 %s %s 1~100" % (cond, ev), d["k1"], d["n1"], cat_q1)
        audit("Q1 %s %s 101~200" % (cond, ev), d["k2"], d["n2"], cat_q1)
r2 = json.load(open("_RESULT_rewind2.json"))
for cond in ("space", "energy", "mat"):
    c = r2["conds"][cond]
    n = c["n"]
    audit("되감기2 Q2′ %s 끝 우세 공유(%s)" % (cond, c["q2_id"]["modal"]), c["q2_id"]["k"], n, cat_q2)
    audit("되감기2 Q3′ %s 가장 흔한 길" % cond, c["q3_path"]["k"], n, cat_q3)
try:
    r8 = json.load(open("_RESULT_rewind8.json"))
    audit("해부8 A2 재료 먼저 끝 우세 공유(rascld)", r8["A2"]["ff"]["top"][0][1], r8["n"]["ff"], cat_q2)
    audit("해부8 A1 덧붙임 racld 몫 < 0.5 (상한)", round(r8["A1"]["ff"] * r8["n"]["ff"]), r8["n"]["ff"], lambda k, n: "흔치 않음" if wilson(k, n)[1] < 0.5 else "(상한 ≥ 0.5)")
except Exception as e:
    print("해부 8 항목 건너뜀:", e)

print("== 문턱 민감도 감사 (🟡 사후 · 판정 아님) — Wilson 범주 판정 · 접시 하나로 뒤집히나")
print("   %-44s %7s  %-14s  %-11s  %-11s  %-11s  %s" % ("판정", "k/n", "Wilson", "지금", "k−1", "k+1", "±1 뒤집힘"))
nflip = 0
for label, k, n, lo, hi, c, cm, cp, flip in rows:
    nflip += flip
    print("   %-44s %7s  [%.3f, %.3f]  %-11s  %-11s  %-11s  %s" % (label, "%d/%d" % (k, n), lo, hi, c, cm, cp, "🚩" if flip else ""))
print("\n   범주 판정 %d개 중 접시 하나로 뒤집히는 것 %d개" % (len(rows), nflip))
print("\n== 부트스트랩 · 비율 문턱 판정 — 구간 끝과 문턱의 거리 (노트에 적힌 값 · 손으로 옮김)")
for label, est, lo, hi, thr, rule in [
    ("해부7 F4 줄어든 몫 R ≥ 1/2", 0.91, 0.87, 0.95, 0.5, "하한 − 문턱"),
    ("해부8 B1 꾸준한 몫 S ≥ 2/3", 1.11, 1.02, 1.22, 2 / 3, "하한 − 문턱"),
    ("해부8 B3a 구성 몫 비 W_age ≥ 2/3", 0.95, 0.75, 1.20, 2 / 3, "하한 − 문턱"),
    ("해부8 K2 ln F / ln I ≥ 2/3 (점추정)", 0.99, float("nan"), float("nan"), 2 / 3, "점추정 − 문턱"),
    ("해부8 F1(해부7) 부풀림 < 1.5 (점추정 최대)", 0.82, float("nan"), float("nan"), 1.5, "문턱 − 점추정"),
    ("해부9 M1 R_mem ≥ 1/2", 0.87, 0.83, 0.91, 0.5, "하한 − 문턱"),
    ("해부5B μ 기울기 하한 > 0", 20.46, 9.91, 32.26, 0.0, "하한 − 문턱"),
    ("해부7 G2 ln 차 상한 < 0", -0.011, -0.021, -0.002, 0.0, "문턱 − 상한"),
    ("해부8 B3b 나이별 몫 하한 > 0", 0.037, 0.006, 0.070, 0.0, "하한 − 문턱"),
    ("해부8 C 개정 뒤 K0 |ln| ≤ 0.10 (최대)", 0.054, float("nan"), float("nan"), 0.10, "문턱 − 최대"),
]:
    if rule.startswith("하한"):
        m = lo - thr
    elif rule.startswith("문턱 − 상한"):
        m = thr - hi
    elif rule.startswith("점추정"):
        m = est - thr
    elif rule.startswith("문턱 − 점추정") or rule.startswith("문턱 − 최대"):
        m = thr - est
    else:
        m = abs(hi)
    print("   %-46s 추정 %+8.3f [%+.3f, %+.3f] · 문턱 %.3f · 여유 %+.3f%s" % (label, est, lo, hi, thr, m, "  🚩 좁음" if abs(m) < 0.05 else ""))
print("\n읽기(🟡): 🚩 는 '문턱을 조금 옮기면 판정이 바뀐다' 는 뜻이지 판정이 틀렸다는 뜻이 아니다. 문턱은 Claude 가 정했다 — 사용자 검토용.")
