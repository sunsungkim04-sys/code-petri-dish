import sys, os, math, random
sys.path.insert(0, os.path.expanduser("~/petri"))
sys.argv = ["x", os.path.expanduser("~/petri")]
import io, contextlib
import pairs25_analyze as m
cap = {}
class Stop(Exception): pass
Orig = m.LSet
class Cap(Orig):
    def __init__(self, name, delta):
        Orig.__init__(self, name, delta); cap[name] = self
        if name == "d8pool": raise Stop()
m.LSet = Cap
buf = io.StringIO()
try:
    with contextlib.redirect_stdout(buf): m.judge()
except Stop: pass
SD, S8 = cap["D"], cap["d8pool"]
def pu_for(key):
    out = []
    for nm, S in (("D", SD), ("d8pool", S8)):
        r = random.Random("%s|ladder|%s" % (key, nm))
        idx = [[r.randrange(S.n) for _ in range(S.n)] for _ in range(m.N_BOOT)]
        out.append([m.flip_cens([S.mean("rascld", mu, ix) for mu in m.MUS]) for ix in idx])
    bD, b8 = out
    diffs = [x - y for x, y in zip(bD, b8) if not (math.isinf(x) and math.isinf(y))]
    return sum(1 for x in diffs if x > 0) / m.N_BOOT, m.pct(diffs, 0.025)
# check frozen stream reproduces
fD, bD = SD.flips(); f8, b8 = S8.flips()
d = [x - y for x, y in zip(bD, b8) if not (math.isinf(x) and math.isinf(y))]
print("frozen P_up %.4f d_lo %+.5f%%p" % (sum(1 for x in d if x > 0) / m.N_BOOT, 100*m.pct(d, 0.025)))
res = [pu_for("haebu25alt%d" % k) for k in range(200)]
ps = sorted(p for p, _ in res)
print("alt200 P_up min %.4f med %.4f max %.4f  >=0.975: %d/200  rank of frozen(<=): %d" % (ps[0], ps[100], ps[-1], sum(p >= 0.975 for p in ps), sum(p <= 0.9855 for p in ps)))
print("alt d_lo min %+.4f med %+.4f max %+.4f (%%p)" % tuple(100*x for x in (min(l for _, l in res), sorted(l for _, l in res)[100], max(l for _, l in res))))
