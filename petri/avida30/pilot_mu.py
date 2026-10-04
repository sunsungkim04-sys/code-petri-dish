from multiprocessing import Pool
import a30_e0 as E
E.ROOT="/home1/minseo1101/projects/avida-port/e30/pilot/mu"
jobs=[(L,200,rn,mu,9001,10000) for L in (3,6) for rn in ("df","ff","do") for mu in (0.0003,0.001)]
with Pool(12) as P:
    for r in P.imap(E.job_de,jobs):
        print(r["L"],r["rule"],r["mu"],"A %.3f R %.3f TR %.3f W %.2f F %.3f supply %.4f I %.3f alive %d viol %s"%(r["A"],r["R"],r["TR"],r["W"],r["F"],r["supply"],r["I"],r["alive_mean"],r["viol"]),flush=True)
