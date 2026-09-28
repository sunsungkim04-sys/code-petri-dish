# -*- coding: utf-8 -*-
"""해부 6 설계용 탐색 요약 — 판정 아님(시드 591~593)."""
import glob, json, statistics as st
print("== A 계보 (μ · 코드별 시드 3개 중앙값)")
print("  %-7s %6s | %6s %6s %6s %6s | %6s %6s | %7s %7s | %s" % ("코드", "μ", "kWT", "k1", "k2", "k3+", "r1", "s1", "WT몫끝", "d≥2몫", "변종출생/WT부모 · 되돌아옴"))
for code in ("racld", "rascld"):
    for mu in ("0p001", "0p005", "0p01"):
        zs = [json.load(open(f)) for f in sorted(glob.glob("pilot6/lin/lin_%s_mu%s_s*.json" % (code, mu)))]
        def k(z, c):
            n, s, _ = z["classes"][c]; return s / n if n else float("nan")
        kw = st.median(k(z, "wt") for z in zs); k1 = st.median(k(z, "d1") for z in zs)
        k2 = st.median(k(z, "d2") for z in zs); k3 = st.median(k(z, "d3") for z in zs)
        s1 = st.median(z["classes"]["d1"][2] / max(z["classes"]["d1"][0], 1) for z in zs)
        last = [z["series"][-1] for z in zs]
        wtf = st.median(r[2] / r[1] for r in last)
        deep = st.median((r[4] + r[5]) / max(r[3] + r[4] + r[5], 1) for r in last)
        u = st.median(z["wt_parent_mut_births"] / z["wt_parent_births"] for z in zs)
        back = st.median(z["wt_from_mut"] / max(z["wt_births"], 1) for z in zs)
        n1 = st.median(z["classes"]["d1"][0] for z in zs)
        print("  %-7s %6s | %6.3f %6.3f %6.3f %6.3f | %6.3f %6.3f | %6.1f%% %6.1f%% | %.3f · %.4f  (d1 n=%d, 놓침 %d)"
              % (code, mu, kw, k1, k2, k3, k1 / kw, s1, 100 * wtf, 100 * deep, u, back, n1, max(z["missed_births"] for z in zs)))

print("\n== B 시간 예산 (시드 3개 합)")
def agg(files, sel):
    tot = {}
    for f in files:
        z = json.load(open(f))
        popend = z["pop_end"]
        for b in z["bins"]:
            if sel(b, popend):
                for kk, v in b.items():
                    if isinstance(v, (int, float)) and kk not in ("t0", "pop0"):
                        tot[kk] = tot.get(kk, 0) + v
    return tot
for tag in ["mat_d04", "mat_d08", "mat_d14", "mat_d32", "space"]:
    for code in ("racld", "rascld", "rsacld"):
        fs = sorted(glob.glob("pilot6/bud/bud_%s_%s_s*.json" % (tag, code)))
        if not fs: continue
        z0 = json.load(open(fs[0]))
        grow_bins = [b["t0"] for b in z0["bins"] if b["pop0"] < 0.5 * z0["pop_end"]]
        for name, sel in (("채움(개체<끝의 50%)", lambda b, pe: b["pop0"] < 0.5 * pe), ("다 찬 뒤(1.5만~)", lambda b, pe: b["t0"] >= 15000)):
            t = agg(fs, sel)
            if not t.get("births"): 
                print("  %-8s %-7s %-18s 출생 0" % (tag, code, name)); continue
            B = t["births"]; S = t["stepped"]
            print("  %-8s %-7s %-18s 틱/출생 %6.1f | 헛손질/출생 %6.1f · 막힘/출생 %6.1f · 복사성공/출생 %5.2f · 헛손질률 %.3f · 막힘률 %.3f | 헛손질몫 %.2f 막힘몫 %.2f"
                  % (tag, code, name, S / B, t["copy_stall"] / B, t["alloc_blocked"] / B, t["copy_ok"] / B,
                     t["copy_stall"] / max(t["copy_stall"] + t["copy_ok"], 1), t["alloc_blocked"] / max(t["alloc_blocked"] + t["alloc_ok"], 1),
                     t["copy_stall"] / S, t["alloc_blocked"] / S))
        print("  %-8s %-7s   (채움 칸 %d개: %s…)" % (tag, code, len(grow_bins), grow_bins[:4]))
