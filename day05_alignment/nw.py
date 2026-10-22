# Needleman-Wunsch score table (global alignment), scores only for now.

MATCH, MISMATCH, GAP = 1, -1, -3   # same scoring as our hand example


def score_table(seq1, seq2):
    """Return the filled score table: seq1 = rows, seq2 = columns."""
    n, m = len(seq1), len(seq2)

    # Step 1: (n+1) x (m+1) zeros; the +1 is the "empty prefix" row/column
    table = [[0] * (m + 1) for _ in range(n + 1)]

    # Step 2: first column and first row hold only gaps (0, -2, -4, ...)
    for i in range(1, n + 1):
        table[i][0] = i * GAP
    for j in range(1, m + 1):
        table[0][j] = j * GAP

    # Step 3 (YOU): fill every other cell
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i-1] == seq2[j - 1]:
                pair = MATCH
            else:
                pair = MISMATCH
            diag = table[i - 1][j - 1] + pair
            up   = table[i - 1][j] + GAP
            left = table[i][j - 1] + GAP

            table[i][j] = max(diag, up, left)      
            

    return table

def trace_back(table, seq1, seq2):
    """Walk from the bottom-right cell to the top-left; return two aligned strings."""
    i, j = len(seq1), len(seq2)   # start at the bottom-right cell
    top, bottom = [], []          # symbols collected while walking backwards

    while i > 0 or j > 0:
        if i > 0 and j > 0:
            if seq1[i - 1] == seq2[j - 1]:
                pair = MATCH
            else:
                pair = MISMATCH
        if i > 0 and j > 0 and table[i][j] == table[i - 1][j - 1] + pair:
            # diagonal: a letter from each sequence
            top.append(seq1[i - 1])
            bottom.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and table[i][j] == table[i - 1][j] + GAP:
            # up: the row letter against a gap
            top.append(seq1[i - 1])
            bottom.append("-")
            i -= 1
        else:
            # left: a gap against the column letter
            top.append("-")
            bottom.append(seq2[j - 1])
            j -= 1    

    return "".join(reversed(top)), "".join(reversed(bottom))    


if __name__ == "__main__":
    s1, s2 = "ACT", "AGT"
    for row in score_table("ACT", "AGT"):
        print(row)
    top, bottom = trace_back(score_table(s1, s2), s1, s2)
    print(top)
    print(bottom)