import numpy as np
import matplotlib
matplotlib.use("Agg")                 
import matplotlib.pyplot as plt

rng = np.random.default_rng(42)        # seed: same numbers every run
n_genes, n_reps = 2000, 20            
phi = 0.1                              

# true mean per gene, evenly spread on a log scale from 10 to 10,000
true_mu = 10 ** rng.uniform(1, 4, size=n_genes)

# TODO 1: Poisson counts, shape (n_genes, n_reps)
pois = rng.poisson(lam=true_mu[:, None], size=(n_genes, n_reps))

# TODO 2: negative binomial counts, same shape
size_n = 1 / phi
p = size_n / (size_n + true_mu[:, None])
nb = rng.negative_binomial(size_n, p, (n_genes, n_reps))

# TODO 3: per-gene mean and sample variance ACROSS replicates, for both matrices
mean_p,  var_p  = pois.mean(axis=1), pois.var(axis=1, ddof=1)
mean_nb, var_nb = nb.mean(axis=1),   nb.var(axis=1, ddof=1)
print(pois.shape)
print(mean_p.shape, mean_p.min(), mean_p.max())
print(mean_nb.shape, mean_nb.min(), mean_nb.max())

# ---- plotting (given) ----
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharex=True, sharey=True)
for ax, m, v, title in [(axes[0], mean_p, var_p, "Poisson"),
                        (axes[1], mean_nb, var_nb, "Negative binomial, phi=0.1")]:
    ok = (m > 0) & (v > 0)             # log axes cannot show zeros
    ax.scatter(m[ok], v[ok], s=4, alpha=0.4)
    ax.plot([10, 1e4], [10, 1e4], color="red", label="variance = mean")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Mean count"); ax.set_title(title)
axes[0].set_ylabel("Variance"); axes[0].legend()
plt.tight_layout(); plt.savefig("overdispersion.png", dpi=150)

near100 = (mean_nb > 80) & (mean_nb < 120)
print("median var/mean near 100, NB:", np.median(var_nb[near100] / mean_nb[near100]))
near100p = (mean_p > 80) & (mean_p < 120)
print("median var/mean near 100, Poisson:", np.median(var_p[near100p] / mean_p[near100p]))