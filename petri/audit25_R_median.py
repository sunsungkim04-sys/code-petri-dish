import sys, os, io, contextlib
sys.path.insert(0, os.path.expanduser("~/petri")); sys.argv = ["x", os.path.expanduser("~/petri")]
import pairs25_analyze as m
orig = m.ib_summary; got = []
class Stop(Exception): pass
def w(IB, codes, name):
    r = orig(IB, codes, name); got.append((name, r[0], r[1], r[3]))
    if len(got) == 2: raise Stop()
    return r
m.ib_summary = w
try:
    with contextlib.redirect_stdout(io.StringIO()): m.judge()
except Stop: pass
import statistics
for name, med, pooled, ns in got:
    v = sorted(pooled.values())
    print(name, "n_bg", ns, "upper-median %.3f" % med, "true median %.3f" % statistics.median(v), " ".join("%s:%.3f" % kv for kv in sorted(pooled.items(), key=lambda x: x[1])))
