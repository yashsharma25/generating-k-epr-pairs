import numpy as np
import itertools
import networkx as nx

import pandas as pd
import re

# example of "4-pairable" graph state:
# one can use LOCC to create any pattern of two EPR states on any subset of 4 qubits

# number of qubits

n = 9

# Argument: binary matrix A
# Returns a matrix whose rows span the nullspace of A over the binary field GF(2)

def null2(A):
    rows, n = A.shape
    X = np.identity(n, dtype=int)
    for i in range(rows):
        y = np.dot(A[i, :], X) % 2
        not_y = (y + 1) % 2
        good = X[:, np.nonzero(not_y)]
        good = good[:, 0, :]
        bad = X[:, np.nonzero(y)]
        bad = bad[:, 0, :]
        if bad.shape[1] > 0:
            bad = np.add(bad,  np.roll(bad, 1, axis=1)) % 2
            bad = np.delete(bad, 0, axis=1)
            X = np.concatenate((good, bad), axis=1)
    return np.transpose(X)


def rank2(A):
    KerA = null2(A)
    return A.shape[1] - KerA.shape[0]

# Argument: stabilizer tableaux of size s x 2n
# each row is a stabilizer, columns range(n) = X part, columns range(n,2n) = Z part
# checks whether qubits a and b are maximally entangled with each other


def test_epr(G, a, b):
    Sa = rank2(G[:, [a, a+n]])
    Sb = rank2(G[:, [b, b+n]])
    Sab = rank2(G[:, [a, b, a+n, b+n]])
    if Sa == 2 and Sb == 2 and Sab == 2:
        return True
    else:
        return False

####################################################################

# define adjacency matrix of the graph state


def get_adjacency_matrix():
    A = np.zeros((n, n), dtype=int)

    for i in range(n):
        j = (i+1) % n
        A[i, j] = 1
        A[j, i] = 1
        j = (i+int(n/2)) % n
        A[i, j] = 1
        A[j, i] = 1

    return A

# compute all pairings of the set {0,1,...,n-1}


def get_n_elem_pairings():
    S4 = itertools.permutations(range(4))
    four_elem_pairings = [sigma for sigma in S4 if sigma[0] <
        sigma[1] and sigma[2] < sigma[3] and sigma[0] < sigma[2]]
    n_elem_pairings = []
    for c in itertools.combinations(range(n), 4):
        for p in four_elem_pairings:
            n_elem_pairings.append([(c[p[0]], c[p[1]]), (c[p[2]], c[p[3]])])
    assert(len(n_elem_pairings) == int(3*n*(n-1)*(n-2)*(n-3)/24))
    return n_elem_pairings

# compute a list of all Pauli measurement bases


def get_pauli_bases():
    xyz = []
    for q in range(n-4):
        xyz.append(['Z', 'X', 'Y'])
    pauli_basis = []
    for elem in itertools.product(*xyz):
        pauli_basis.append(elem)

    # pauli_basis now contains all possible measurements of the 6 qubits
    return pauli_basis


def get_stabilizer_tableaux(A):
    # stabilizer tableaux of the graph state
    G = np.zeros((n, 2*n), dtype=int)
    for i in range(n):
        G[i, i] = 1
        G[i, n:2*n] = A[i, :]
    GT = np.transpose(G)
    return G, GT

# symplectric inner product matrix
def get_symplectric_inner_product():
    Lambda = np.zeros((2*n, 2*n), dtype=int)
    for q in range(n):
        Lambda[q][q+n] = 1
        Lambda[q+n][q] = 1
    return Lambda


def epr_generation_checker(A):
    n_elem_pairings = get_n_elem_pairings()
    pauli_basis = get_pauli_bases()
    generated_pairs = 0

    # A = get_adjacency_matrix()
    G, GT = get_stabilizer_tableaux(A)

    Lambda = get_symplectric_inner_product()

    for pp in n_elem_pairings:
        isOK = False
        for meas in pauli_basis:
            a, b = pp[0]
            c, d = pp[1]
            qubits_to_measure = np.sort(list(set(range(n))-set([a, b, c, d])))

            # stabilizer tableaux for the measured qubits
            GM = np.zeros((n-4, 2*n), dtype=int)
            for i in range(n-4):
                q = qubits_to_measure[i]
                if meas[i] == 'X':
                    GM[i, q] = 1
                if meas[i] == 'Y':
                    GM[i, q] = 1
                    GM[i, q+n] = 1
                if meas[i] == 'Z':
                    GM[i, q+n] = 1

            # find graph state stabilizers commuting with all measurements
            Gcom = np.dot(null2(np.dot(GM, np.dot(Lambda, GT)) % 2), G) % 2

            if test_epr(Gcom, a, b) and test_epr(Gcom, c, d):
                print('EPR pairs', pp, 'measurement basis=', end='')
                generated_pairs += 1
                cnt = 0
                for q in range(n):
                    if q in qubits_to_measure:
                        print(meas[cnt], sep='', end='')
                        cnt = cnt + 1
                    else:
                        print('*', sep='', end='')
                print('')
                isOK = True
                break

        if not(isOK):
            print('EPR pairs ', pp, ' cannot be generated')

            # Should we terminate the search for each class
            print('Terminating the search...')
            #exit()
            return generated_pairs

    return generated_pairs

def main():
    df = pd.read_csv('orbits-data/9qubitorbitsCi_labeled.csv')

    # Each row is an entanglement class
    for row in zip(df['i'], df['graph_list']):
        i, graph_list = row

        graphs = graph_list.split('),')

        # remove the leftmost '('
        graphs[0] = graphs[0][1:]

        # remove the two rightmost ')'s
        graphs[-1] = graphs[-1][:-2]

        print("Graphs in class ", i, " = ", len(graphs))
        for index,g in enumerate(graphs):
            g = re.sub(r'\s+', '', g)            
            edges = g[1:].split(',')

            A = np.zeros((n, n), dtype=int)

            for e in edges:
                edge_start = int(e.split('-')[0]) - 1
                edge_end = int(e.split('-')[1]) - 1

                # create adjacency matrix
                A[edge_start][edge_end] = 1
                A[edge_end][edge_start] = 1

            generated_pairs = epr_generation_checker(A)

            if generated_pairs == int(3*n*(n-1)*(n-2)*(n-3)/24):
                print("Generated all Pairs =", generated_pairs, "for clas = ", i, " and graph = ", index)
                print("G =", g)
                print(A)


            # print(A)
            # print("____________________________________")

        # print(
        #     "******************************************************************************")
    
        

def single_class():
    df = pd.read_csv('orbits-data/9qubitorbitsCi_labeled.csv')

    # Each row is an entanglement class
    for row in zip(df['i'], df['graph_list']):
        i, graph_list = row

        graphs = graph_list.split('),')

        # remove the leftmost '('
        graphs[0] = graphs[0][1:]

        # remove the two rightmost ')'s
        graphs[-1] = graphs[-1][:-2]

        print("Graphs in class ", i, " = ", len(graphs))
        g = graphs[0]
        g = re.sub(r'\s+', '', g)            
        edges = g[1:].split(',')

        A = np.zeros((n, n), dtype=int)

        for e in edges:
            edge_start = int(e.split('-')[0]) - 1
            edge_end = int(e.split('-')[1]) - 1

            # create adjacency matrix
            A[edge_start][edge_end] = 1
            A[edge_end][edge_start] = 1


        generated_pairs = epr_generation_checker(A)
        print("Generated Pairs =", generated_pairs)

        if generated_pairs == int(3*n*(n-1)*(n-2)*(n-3)/24):
            print("Generated all Pairs =", generated_pairs, "for class = ", i, " and graph = ", index)
            print("G =", g)
            print(A)

        # print(
        #     "******************************************************************************")
    
    print("Total Generated Pairs =", generated_pairs)

#main()

single_class()