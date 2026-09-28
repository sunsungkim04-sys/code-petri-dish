#!/usr/bin/env python3
"""해부 11 계측기 패치 (09-20). budget.js v1.1.0 → v1.2.0 · dms.js 에 --age0/--age-var.
기본값이면 궤적 · 출력(더한 필드 제외)이 같아야 한다 — reg11_check.py 가 대조한다. 한 번만 돌린다."""
import io, sys

def patch(path, pairs):
    t = io.open(path, encoding="utf-8").read()
    for a, b in pairs:
        if t.count(a) != 1:
            sys.exit("패치 자리 %d곳: %s · %r" % (t.count(a), path, a[:60]))
        t = t.replace(a, b)
    io.open(path, "w", encoding="utf-8").write(t)
    print("OK", path)

patch("budget.js", [
    ("const BUDGET_VERSION = '1.1.0';", "const BUDGET_VERSION = '1.2.0';"),
    (" * 사용: node budget.js --cond mat",
     " * v1.2.0(해부 11): --age0 N --age-var N(sim 0.3.5 · 세계 전체) · 칸마다 죽음 수 · 죽은 나이 합 · 칸 시작의 글자별 자유 재료를 더 적는다.\n"
     " *         읽기만 한다 — 옵션을 안 주면 궤적 · 옛 필드가 v1.1.0 과 같다(reg11_check.py).\n"
     " *\n * 사용: node budget.js --cond mat"),
    ("const outDir = args.out || '.';\nfor (const ch of code)",
     "const outDir = args.out || '.';\n"
     "const age0V = args.age0 ? parseInt(args.age0, 10) : null;\n"
     "const ageVarV = args['age-var'] ? parseInt(args['age-var'], 10) : null;\n"
     "if ((age0V !== null || ageVarV !== null) && Sim.VERSION !== '0.3.5') { console.error('--age0 은 sim 0.3.5 가 필요하다'); process.exit(2); }\n"
     "for (const ch of code)"),
    ("const w = Sim.makeWorld(Object.assign({ seed, density, mu: 0, light: 1, ancestor: code }, COND[cond]));",
     "const wopts = Object.assign({ seed, density, mu: 0, light: 1, ancestor: code }, COND[cond]);\n"
     "if (age0V !== null) wopts.age0 = age0V;\nif (ageVarV !== null) wopts.ageVar = ageVarV;\n"
     "const w = Sim.makeWorld(wopts);"),
    ("const mk = (t, pop) => { const o = { t0: t, pop0: pop };",
     "const freeNow = () => { if (!w.mat) return null; const f = new Array(Sim.K).fill(0); const m = w.mat; for (let j = 0; j < m.length; j++) f[j % Sim.K] += m[j]; return f; };\n"
     "const mk = (t, pop) => { const o = { t0: t, pop0: pop, deaths: 0, death_age_sum: 0, free0: freeNow() };"),
    ("  W.births += w.births - b0;\n",
     "  W.births += w.births - b0;\n"
     "  for (let q = 0; q < n; q++) if (pre[q].dead) { W.deaths++; W.death_age_sum += pre[q].age; }\n"),
    ("  budget_version: BUDGET_VERSION, sim_version: Sim.VERSION, cond,",
     "  budget_version: BUDGET_VERSION, sim_version: Sim.VERSION, opts: wopts, mat_total0: w.matTotal0, cond,"),
    ("const tag = `bud_${cond}${COND[cond].mat ? '_d' + String(density).padStart(2, '0') : ''}_${code}_s${String(seed).padStart(5, '0')}`;",
     "const tag = `bud_${cond}${COND[cond].mat ? '_d' + String(density).padStart(2, '0') : ''}_${code}_s${String(seed).padStart(5, '0')}${age0V !== null ? '_a' + age0V : ''}`;"),
])

patch("dms.js", [
    ("const ancestor = args.ancestor || null;\n",
     "const ancestor = args.ancestor || null;\n"
     "// 해부 11(09-20 추가): --age0 N --age-var N 이면 **배경과 팔 모두** 그 수명으로 돈다(sim 0.3.5 · 세계 전체 옵션). 기본값(안 줌)이면 출력이 전과 같다.\n"
     "const age0V = args.age0 ? parseInt(args.age0, 10) : null;\n"
     "const ageVarV = args['age-var'] ? parseInt(args['age-var'], 10) : null;\n"
     "if ((age0V !== null || ageVarV !== null) && Sim.VERSION !== '0.3.5') { console.error('--age0 은 sim 0.3.5 가 필요하다'); process.exit(2); }\n"),
    ("if (ancestor) opts.ancestor = ancestor;\n",
     "if (ancestor) opts.ancestor = ancestor;\nif (age0V !== null) opts.age0 = age0V;\nif (ageVarV !== null) opts.ageVar = ageVarV;\n"),
    ("if (ancestor) out.ancestor = ancestor;\n",
     "if (ancestor) out.ancestor = ancestor;\nif (age0V !== null) out.age0 = age0V;\nif (ageVarV !== null) out.age_var = ageVarV;\n"),
])
