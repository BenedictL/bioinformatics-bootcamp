"""Day 3: minimal FASTQ parser (no Biopython).

Run in WSL with:
    XDG_RUNTIME_DIR=/mnt/wslg/runtime-dir python parse_fastq.py
"""
import matplotlib.pyplot as plt


def phred_scores(qual):
    """Return a list of Phred Q values, one per character."""
    scores = []
    for c in qual:
        q = ord(c) - 33
        scores.append(q)
    return scores


def gc_fraction(seq):
    """Return the fraction of G and C bases in one sequence."""
    count = 0
    for base in seq:
        if base == "G" or base == "C":
            count += 1
    return count / len(seq)


totals = [0] * 10
counts = [0] * 10

with open("toy.fastq") as f:
    while True:
        header = f.readline().strip()
        if not header:
            break
        seq = f.readline().strip()
        plus = f.readline().strip()
        qual = f.readline().strip()
        # print(header, gc_fraction(seq), len(seq))
        scores = phred_scores(qual)
        for i, q in enumerate(scores):
            totals[i] += q
            counts[i] += 1

means = []
for i in range(10):
    means.append(round(totals[i] / counts[i], 2))
print(means)

x = range(1, 11)
y = means

plt.plot(x, y, marker="o")
plt.axhline(30, linestyle="--")
plt.xlabel("position in read")
plt.ylabel("mean Phred Q")
plt.title("Mean quality per position")
plt.savefig("mean_quality.png", dpi=150)
plt.show()
