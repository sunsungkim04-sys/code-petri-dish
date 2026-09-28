#!/usr/bin/env python3
"""해부 12 계측기 패치 (09-20). dms.js 에 --kids 1(팔 안에서 표지 계통의 출생을 부모 코드 · 자식 코드로 센다 · 읽기만) 과
--arms nbr --nbr-of CODE(임의 코드의 한 글자 이웃 전부 + 그 코드 자신을 WT 공동체에 끼운다)를 더한다.
기본값이면 출력이 전과 같아야 한다 — reg12_check.py 가 대조한다. 한 번만 돌린다."""
import io, sys

def patch(path, pairs):
    t = io.open(path, encoding="utf-8").read()
    for a, b in pairs:
        if t.count(a) != 1:
            sys.exit("패치 자리 %d곳: %s · %r" % (t.count(a), path, a[:60]))
        t = t.replace(a, b)
    io.open(path, "w", encoding="utf-8").write(t)
    print("OK", path)

patch("dms.js", [
    ("const noMat = args['no-mat'] === '1';\n",
     "const noMat = args['no-mat'] === '1';\n"
     "// 해부 12(09-20 추가): --kids 1 이면 팔을 돌리는 동안 **표지 계통의 출생을 밖에서 읽어** 부모가 끼운 코드 그대로일 때의 자식 코드 분포를 적는다\n"
     "// (틱 전에 '자식 칸 → 부모 코드' 를 적어 두고 틱 뒤에 새 id 를 그 칸에서 찾는다 · 읽기만 — 궤적 · 옛 필드 동일).\n"
     "// 복사 없이 코드가 바뀐 것(우주선)도 따로 센다. --arms nbr --nbr-of CODE 면 CODE 의 한 글자 이웃 전부와 CODE 자신을 끼운다.\n"
     "const kidsOn = args.kids === '1';\n"
     "const nbrOf = args['nbr-of'] || null;\n"
     "if (armsMode === 'nbr' && !nbrOf) { console.error('--arms nbr 에는 --nbr-of CODE 가 필요하다'); process.exit(2); }\n"),
    ("  const series = [[0, ...count()]];\n  let stopped = -1;\n  for (let t = 1; t <= T; t++) {\n    Sim.tick(w);\n",
     "  const series = [[0, ...count()]];\n  let stopped = -1;\n"
     "  const KD = kidsOn ? { b_exact: 0, b_other: 0, unknown: 0, kids: new Map(), cosmic: new Map() } : null;\n"
     "  const bump = (m, k) => m.set(k, (m.get(k) || 0) + 1);\n"
     "  for (let t = 1; t <= T; t++) {\n"
     "    let pre = null, id0 = 0, exactList = null;\n"
     "    if (kidsOn) {\n"
     "      pre = new Map(); id0 = w.nextId; exactList = [];\n"
     "      for (const o of w.alive) if (o.tag === 1) { if (o.child) pre.set(o.child.i, o.key); if (o.key === key) exactList.push(o); }\n"
     "    }\n"
     "    Sim.tick(w);\n"
     "    if (kidsOn) {\n"
     "      for (const o of w.alive) if (o.id >= id0 && o.tag === 1) {\n"
     "        const pk = pre.get(o.i);\n"
     "        if (pk === undefined) KD.unknown++;\n"
     "        else if (pk === key) { KD.b_exact++; if (o.key !== key) bump(KD.kids, o.key); }\n"
     "        else KD.b_other++;\n"
     "      }\n"
     "      for (const o of exactList) if (!o.dead && o.key !== key) bump(KD.cosmic, o.key);\n"
     "    }\n"),
    ("  return { kind, key, labels, post_seed: postSeed, inj, cons_ok: consAfterInj && consEnd, stopped, series, elapsed_ms: Date.now() - ta };",
     "  const ret = { kind, key, labels, post_seed: postSeed, inj, cons_ok: consAfterInj && consEnd, stopped, series, elapsed_ms: Date.now() - ta };\n"
     "  if (kidsOn) ret.kids = { b_exact: KD.b_exact, b_other: KD.b_other, unknown: KD.unknown, kids: [...KD.kids.entries()].sort((a, b) => b[1] - a[1]), cosmic: [...KD.cosmic.entries()].sort((a, b) => b[1] - a[1]) };\n"
     "  return ret;"),
    ("    : armsMode === 'list' ? muts.filter(([key]) => codes.includes(key))\n",
     "    : armsMode === 'list' ? muts.filter(([key]) => codes.includes(key))\n"
     "    : armsMode === 'nbr' ? [[nbrOf, ['self']], ...mutantList(nbrOf).filter(([key]) => key !== WT)]\n"),
    ("if (ancestor) out.ancestor = ancestor;\n",
     "if (ancestor) out.ancestor = ancestor;\nif (kidsOn) out.kids_on = true;\nif (nbrOf) out.nbr_of = nbrOf;\n"),
    ("const tag = `dms_${cond}_${WT}_s${String(seed).padStart(5, '0')}_${armsMode}`;",
     "const tag = `dms_${cond}_${WT}_s${String(seed).padStart(5, '0')}_${armsMode}${nbrOf ? '_' + nbrOf : ''}`;"),
])
