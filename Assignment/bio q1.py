def global_alignment_affine_gap(v, w, score_matrix, gap_opening, gap_extension):
    def score(x, y):
        return score_matrix[x][y]

    n, m = len(v), len(w)
    M = [[0] * (m + 1) for _ in range(n + 1)]
    I_x = [[float('-inf')] * (m + 1) for _ in range(n + 1)]
    I_y = [[float('-inf')] * (m + 1) for _ in range(n + 1)]
    traceback_matrix = [[None] * (m + 1) for _ in range(n + 1)]

    # Initialization
    for i in range(1, n + 1):
        I_y[i][0] = gap_opening + (i - 1) * gap_extension
        M[i][0] = I_y[i][0]
        traceback_matrix[i][0] = '↑'  # Gap in w
    for j in range(1, m + 1):
        I_x[0][j] = gap_opening + (j - 1) * gap_extension
        M[0][j] = I_x[0][j]
        traceback_matrix[0][j] = '←'  # Gap in v

    # Fill DP tables
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            M[i][j] = max(M[i-1][j-1], I_x[i-1][j-1], I_y[i-1][j-1]) + score(v[i-1], w[j-1])
            I_x[i][j] = max(I_x[i][j-1] + gap_extension, M[i][j-1] + gap_opening)
            I_y[i][j] = max(I_y[i-1][j] + gap_extension, M[i-1][j] + gap_opening)

            if M[i][j] >= I_x[i][j] and M[i][j] >= I_y[i][j]:
                traceback_matrix[i][j] = '↖'  # Match/mismatch
            elif I_x[i][j] > I_y[i][j]:
                traceback_matrix[i][j] = '←'  # Gap in v
            else:
                traceback_matrix[i][j] = '↑'  # Gap in w

    # Traceback
    alignment_v, alignment_w = '', ''
    i, j = n, m
    while i > 0 or j > 0:
        if traceback_matrix[i][j] == '↖':
            alignment_v = v[i-1] + alignment_v
            alignment_w = w[j-1] + alignment_w
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == '←':
            alignment_v = '-' + alignment_v
            alignment_w = w[j-1] + alignment_w
            j -= 1
        elif traceback_matrix[i][j] == '↑':
            alignment_v = v[i-1] + alignment_v
            alignment_w = '-' + alignment_w
            i -= 1

    return M, I_x, I_y, max(M[n][m], I_x[n][m], I_y[n][m]), alignment_v, alignment_w

# Example usage
v = "ACAGT"
w = "ACGT"
score_matrix = {
    'A': {'A': 2, 'C': -1, 'G': -1, 'T': -1},
    'C': {'A': -1, 'C': 2, 'G': -1, 'T': -1},
    'G': {'A': -1, 'C': -1, 'G': 2, 'T': -1},
    'T': {'A': -1, 'C': -1, 'G': -1, 'T': 2}
}
gap_opening = -2
gap_extension = -1

M, I_x, I_y, result, alignment_v, alignment_w = global_alignment_affine_gap(v, w, score_matrix, gap_opening, gap_extension)

print("Highest-scoring global alignment:", result)
print("Alignment v:", alignment_v)
print("Alignment w:", alignment_w)

print("\nMatrix M:")
for row in M:
    print(row)
print("\nMatrix I_x:")
for row in I_x:
    print(row)
print("\nMatrix I_y:")
for row in I_y:
    print(row)
