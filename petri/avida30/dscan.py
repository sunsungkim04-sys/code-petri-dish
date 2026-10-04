import sys, json
from multiprocessing import Pool
import a30lib as A
T=10000
import os
DENS=[int(x) for x in os.environ["DENS"].split()]
def job(a):
    d,L=a
    p=f"/home1/minseo1101/projects/avida-port/e30/pilot/dscan/d{d}_L{L}"
    ev=[f"u begin InjectSequence {A.FAMILY_S[L]} {A.CENTER} {A.CENTER+1} -1 0", f"u {T} Exit"]
    info=A.run(p, ev, {"MATERIAL_DENSITY":d, "RANDOM_SEED":9001, "COPY_DRAW_RULE":0})
    rows=A.read_census(p); s=A.window_sums(rows,0,T//2,T); m=A.measures(s,L,18,0)
    return (d,L,info["elapsed_s"],info["rc"],round(m["A"],3),round(m["TR"],3),round(m["W"],2),round(m["alive_mean"]),m["placed"],A.audit_violations(p))
jobs=[(d,L) for d in DENS for L in (3,6)]
with Pool(8) as P:
    for r in P.imap(job,jobs): print(*r, flush=True)
