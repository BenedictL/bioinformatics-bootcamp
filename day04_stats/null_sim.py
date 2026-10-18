import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")

rng = np.random.default_rng(42)       # seed: same numbers every run

n_genes = 10_000
n_per_group = 3

# TODO 1: make two arrays of shape (n_genes, n_per_group), a and b,
# both drawn from the SAME normal distribution (mean 0, sd 1).
# Hint: rng.normal(...) has a size argument.
a = rng.normal(loc=0, scale=1, size=(n_genes, n_per_group))
b = rng.normal(loc=0, scale=1, size=(n_genes, n_per_group))

# t-test on every row (gene) at once; axis=1 means "test along each row"
pvals = stats.ttest_ind(a, b, axis=1).pvalue

# TODO 2: count the genes with p < 0.05 (hint: a comparison gives True/False, and sum() counts the True values)
n_hits = (pvals < 0.05).sum()
print("genes with p < 0.05:", n_hits)

plt.hist(pvals, bins=20)
plt.xlabel("p-value")
plt.ylabel("number of genes")
plt.savefig("pvalue_hist.png")     

from statsmodels.stats.multitest import multipletests

# [0] of the result is True/False per gene: True = flagged after correction
bonf = multipletests(pvals, alpha=0.05, method="bonferroni")[0]
bh   = multipletests(pvals, alpha=0.05, method="fdr_bh")[0]

print(bonf.sum())
print(bh.sum())


# TODO: print how many genes each method flags (same .sum() trick as before)