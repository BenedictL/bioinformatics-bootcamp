import numpy as np

rng = np.random.default_rng(1)
m = 1000
p = rng.uniform(size=m)                       # 1000 genes, all null at first
p[:20] = rng.uniform(0, 0.0005, size=20)      # genes 0-19 have a real effect (tiny p)
order = np.argsort(p)                         # gene indices, smallest p first
ps = p[order]                                 # the p-values in that sorted order

ranks = np.arange(1, m + 1)
threshold = ranks / m * 0.5                  # BH line: (rank / m) * 0.05
passing = np.where(ps <= threshold)[0]        # sorted positions that are under the line
k = passing.max() + 1 if len(passing) > 0 else 0   # largest rank under the line
called = np.zeros(m, dtype=bool)
called[order[:k]] = True   
print("k (largest rank under the line):", k)                   # compare each p to its line value
print("genes called significant:", passing.sum())
print("called among the 20 real genes:", passing[:20].sum())
print("called among the 980 null genes:", passing[20:].sum())

import matplotlib
matplotlib.use("Agg")                         # save to a file, no window (as in your other scripts)
import matplotlib.pyplot as plt

is_real = order < 20                          # True where the sorted gene is one of genes 0-19
plt.figure(figsize=(7, 4.5))
plt.scatter(ranks[~is_real], ps[~is_real], s=8, color="tab:gray", label="null genes")
plt.scatter(ranks[is_real], ps[is_real], s=20, color="tab:green", label="real effects")
plt.plot(ranks, threshold, color="tab:orange", label="BH line (q = 0.05)")
plt.axvline(k, linestyle="--", color="black", label=f"k = {k}")
plt.xscale("log"); plt.yscale("log")
plt.xlabel("rank (smallest p first)"); plt.ylabel("p-value")
plt.legend(); plt.tight_layout(); plt.savefig("bh_plot.png", dpi=150)