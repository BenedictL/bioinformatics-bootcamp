import numpy as np

rng = np.random.default_rng(42)          # seed: same numbers every run

true_effect = 10      # real treatment effect
batch_shift = 50      # Tuesday adds this to every sample
baseline = 100
noise_sd = 5
n_per_cell = 4

# ---- two designs (given) ----
# balanced: 4 samples in each treatment x day cell
treat_bal = np.repeat([1, 1, 0, 0], n_per_cell)   # 1 = treated
batch_bal = np.repeat([0, 1, 0, 1], n_per_cell)   # 1 = Tuesday
# confounded: all treated on Monday, all controls on Tuesday
treat_con = np.repeat([1, 0], 2 * n_per_cell)
batch_con = 1 - treat_con

def simulate(treat, batch):
    noise = rng.normal(0, noise_sd, size=len(treat))
    return baseline + true_effect * treat + batch_shift * batch + noise
    

y_bal = simulate(treat_bal, batch_bal)   
y_con = simulate(treat_con, batch_con)   
print(y_bal.min().round(1), y_bal.max().round(1))
def ols(y, cols):
    X = np.column_stack([np.ones(len(y))] + cols)
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    # ---- given: rank, residual variance, standard errors ----
    rank = np.linalg.matrix_rank(X)
    resid = y - X @ beta
    sigma2 = resid @ resid / (len(y) - rank)
    se = np.sqrt(sigma2 * np.diag(np.linalg.pinv(X.T @ X)))
    return beta, se, rank

for name, y, t, b in [("balanced", y_bal, treat_bal, batch_bal),
                      ("confounded", y_con, treat_con, batch_con)]:
    beta1, se1, rank1 = ols(y, [t])        # y ~ treatment
    beta2, se2, rank2 = ols(y, [t, b])     # y ~ treatment + batch
    print(name)
    print(f"  ~ treatment:         effect {beta1[1]:.2f}  se {se1[1]:.2f}  rank {rank1}")
    print(f"  ~ treatment + batch: effect {beta2[1]:.2f}  se {se2[1]:.2f}  rank {rank2}  batch {beta2[2]:.2f}")
