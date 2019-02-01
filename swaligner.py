#!/usr/bin/env python

import numpy as np

#penalty box
match    = 3
mismatch = -3
gap      = -2

#sequences. seq2 should always be larger than seq1
seq1 = 'TGTTACGG'
seq2 = 'GGTTGACTA'

#scoring matrix size
rows = len(seq1) + 1
cols = len(seq2) + 1

def create_score_matrix(rows, cols):
    '''
    Create a matrix of scores representing trial alignments of the two sequences.
    Sequence alignment can be treated as a graph search problem. This function
    creates a graph (2D matrix) of scores, which are based on trial alignments
    of different base pairs. The path with the highest cummulative score is the
    best alignment.
    '''
    score_matrix = [[0 for col in range(cols)] for row in range(rows)]
    # Fill the scoring matrix.
    max_score = 0
    max_pos   = None    # The row and columbn of the highest score in matrix.
    for i in range(1, rows):
        for j in range(1, cols):
            score = calc_score(score_matrix, i, j)
            if score > max_score:
                max_score = score
                max_pos   = (i, j)
            score_matrix[i][j] = score
    assert max_pos is not None, 'the x, y position with the highest score was not found'
    return score_matrix, max_pos

def calc_score(matrix, x, y):
    #Calculate score for a given x, y position in the scoring matrix.
    #The score is based on the up, left, and upper-left diagnol neighbors
    similarity = match if seq1[x - 1] == seq2[y - 1] else mismatch
    diag_score = matrix[x - 1][y - 1] + similarity
    up_score   = matrix[x - 1][y] + gap
    left_score = matrix[x][y - 1] + gap
    return max(0, diag_score, up_score, left_score)

def print_matrix(matrix):
    '''
    Print the scoring matrix.
    ex:
    0   0   0   0   0   0
    0   2   1   2   1   2
    0   1   1   1   1   1
    0   0   3   2   3   2
    0   2   2   5   4   5
    0   1   4   4   7   6
    '''
    print(np.matrix(matrix).T)

###add your function(s) to find a solution here.

#from the matrix and the start position returns a list of the traversal
def trace(matrix, start_pos):
    result = [start_pos]
    current_pos = start_pos
    start_val = matrix[int(start_pos[0])][int(start_pos[1])]
    while start_val != 0:
        diag_score = matrix[int(current_pos[0]) -1][int(current_pos[1]) - 1]
        left_score = matrix[int(current_pos[0])][int(current_pos[1] - 1)]
        upper_score = matrix[int(current_pos[0]) -1][int(current_pos[1])]
        if diag_score >= left_score:
            if diag_score >= upper_score:
                result.append((int(current_pos[0]) -1, int(current_pos[1]) - 1))
                current_pos = (int(current_pos[0]) -1, int(current_pos[1]) - 1)
        elif left_score >= upper_score:
            if left_score >= diag_score:
                result.append((int(current_pos[0]), int(current_pos[1]) - 1))
                current_pos = (int(current_pos[0]), int(current_pos[1]) - 1)
        elif upper_score >= left_score:
            if upper_score >= diag_score:
                result.append((int(current_pos[0]) -1, int(current_pos[1])))
                current_pos = (int(current_pos[0]) -1, int(current_pos[1]))
        start_val = matrix[int(current_pos[0])][int( current_pos[1])]
    return result

#from the list of the traversal returns a string representation
def trace_string(trace):
    result = ''
    for pos in range(len(trace)):
        if (pos + 1) != len(trace):
            result = result + 'r' + str(trace[pos][1]) + 'c' + str(trace[pos][0]) + ' -> '
        else:
            result = result + 'r' + str(trace[pos][1]) + 'c' + str(trace[pos][0])
    print(result)

#returns a list of the actual matrix values of the traversal
def check(matrix, trace):
    result = []
    for i in trace:
        result.append(matrix[i[0]][i[1]])
    return result

###end of your function(s)

if __name__ == '__main__':
    #my main
    score_matrix, start_pos = create_score_matrix(rows, cols)
    print_matrix(score_matrix)
    traversal = trace(score_matrix, start_pos)
    trace_string(traversal)
    #print(check(score_matrix, traversal))
