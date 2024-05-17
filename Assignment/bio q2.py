import numpy as np
import matplotlib.pyplot as plt

def init_matrices(seq1, seq2, gap_penalty):
    m, n = len(seq1), len(seq2)
    score_matrix = np.zeros((m+1, n+1))
    traceback_matrix = np.zeros((m+1, n+1), dtype=str)
    
    for i in range(1, m+1):
        score_matrix[i, 0] = score_matrix[i-1, 0] + gap_penalty
        traceback_matrix[i, 0] = '↑'
    
    for j in range(1, n+1):
        score_matrix[0, j] = score_matrix[0, j-1] + gap_penalty
        traceback_matrix[0, j] = '←'
    
    return score_matrix, traceback_matrix

def needleman_wunsch(seq1, seq2, match_score, mismatch_penalty, gap_penalty):
    m, n = len(seq1), len(seq2)
    score_matrix, traceback_matrix = init_matrices(seq1, seq2, gap_penalty)
    
    for i in range(1, m+1):
        for j in range(1, n+1):
            match = score_matrix[i-1, j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            delete = score_matrix[i-1, j] + gap_penalty
            insert = score_matrix[i, j-1] + gap_penalty
            score_matrix[i, j] = max(match, delete, insert)
            
            if score_matrix[i, j] == match:
                traceback_matrix[i, j] = '↖'
            elif score_matrix[i, j] == delete:
                traceback_matrix[i, j] = '↑'
            else:
                traceback_matrix[i, j] = '←'
            
            # Visualization after each step
            visualize_matrix(score_matrix, traceback_matrix, i, j, seq1, seq2)
    
    alignment_a, alignment_b = traceback(seq1, seq2, traceback_matrix)
    return score_matrix, traceback_matrix, alignment_a, alignment_b

def visualize_matrix(score_matrix, traceback_matrix, i, j, seq1, seq2):
    fig, ax = plt.subplots()
    cax = ax.matshow(score_matrix, cmap='viridis')
    plt.title(f"Update at ({i}, {j})")
    fig.colorbar(cax)
    
    # Annotate the matrix with the values
    for x in range(score_matrix.shape[0]):
        for y in range(score_matrix.shape[1]):
            ax.text(y, x, f'{int(score_matrix[x, y])}\n{traceback_matrix[x, y]}', va='center', ha='center', color='red')

    plt.xlabel('Seq2: ' + ', '.join(seq2))
    plt.ylabel('Seq1: ' + ', '.join(seq1))
    plt.show()

def traceback(seq1, seq2, traceback_matrix):
    alignment_a = ""
    alignment_b = ""
    i, j = len(seq1), len(seq2)
    
    while i > 0 or j > 0:
        if traceback_matrix[i, j] == '↖':
            alignment_a = seq1[i-1] + alignment_a
            alignment_b = seq2[j-1] + alignment_b
            i -= 1
            j -= 1
        elif traceback_matrix[i, j] == '↑':
            alignment_a = seq1[i-1] + alignment_a
            alignment_b = '-' + alignment_b
            i -= 1
        else:  # traceback_matrix[i, j] == '←'
            alignment_a = '-' + alignment_a
            alignment_b = seq2[j-1] + alignment_b
            j -= 1
    
    return alignment_a, alignment_b

# Example DNA sequences
seq1 = "GATTACA"
seq2 = "GCATGCU"

# Perform Needleman-Wunsch and visualize each step
score_matrix, traceback_matrix, alignment_a, alignment_b = needleman_wunsch(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-1)

# Print the final alignment
print("Alignment:")
print(alignment_a)
print(alignment_b)
