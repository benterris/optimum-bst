"""
This program computes an optimum binary search tree, implementing
the idea from http://www.inrg.csie.ntu.edu.tw/algorithm2014/presentation/Knuth71.pdf

Here, we included in this implementation the `alpha` part of the article,
which means we take in account the frequencies of words appearing between selected words.

This can be used without the `alpha` frequencies (providing only the `beta` ones to the main function),
or optionally we can provide a list with `alpha` frequencies. In this case the length of the `alpha` list
must be greater than the `beta` list by one unit.
"""

import pprint
import numpy as np


def main(betaList, alphaList=None):
    """
    From a list of frequencies, find the related optimal tree
    We expect a list of frequencies for the words (a list of beta values),
    but one can also optionally provide a list of alpha values (see above).

    :param beta: a list of frequencies for all the words to be put in the tree
    :param alpha: optional list of frequencies for all the words inbetween the latter ones 
    """

    alpha = np.zeros(
        len(betaList) + 1) if alphaList == None else np.array(alphaList)
    beta = np.array(betaList)

    R = findBestRoots(alpha, beta)

    return parentsFromRoots(R)


def computeWeight(alpha, beta):
    """
    Compute W, the total weight matrix

    :param alpha: list of the weight of the alpha parameters
    :param beta: list of the weight of the beta parameters
    """
    W = np.zeros((n, n))
    # Initialise with values for a single node
    for i in range(n):
        W[i, i] = alpha[i] + beta[i] + alpha[i + 1]

    # Compute all the necessary sums, each time adding a new alpha + beta value
    for i in range(n):
        for j in range(i + 1, n):
            W[i, j] = W[i, j - 1] + beta[j] + alpha[j + 1]
    return W


def findBestRoots(alpha, beta):
    """
    Compute R and P, respectively the roots and the weighted paths matrix

    :return R: The roots matrix, containing for R[i, j] the index of the root of the optimal
                subtree built from the words A_i, A_i+1 ... A_j
    """
    # Pre compute all the sums
    W = computeWeight(alpha, beta)

    # Initialisation for P and R
    # We initialise them with -1 to highlight the unused cells
    P = np.ones((n, n)) * -1
    R = (np.ones((n, n)) * -1).astype(int)

    # Compute values for a single node
    for i in range(n):
        # The only candidate to be selected as root is the node itself
        R[i, i] = i
        P[i, i] = alpha[i] + beta[i] + alpha[i + 1]

    # Then for every possible interval length, we compute the P and R value
    # We proceed this way because to compute P and R for a given interval, we
    # only need values computed on shorter intervals
    for groupLen in range(2, n + 1):
        for i in range(n - groupLen + 1):
            # Our aim here is to fill P[i, j] and R[i, j]

            # For convenience, we retrieve j from i and the interval length
            j = i + groupLen - 1

            # We compute all the possible weighted paths for all the possible roots
            # NOTE: we only have to try the roots between R[i, j - 1] and R[i + 1, j], thanks
            # to the findings of the article
            weights = [(P[i, k - 1] if k - 1 >= i else 0) +
                       (P[k + 1, j] if k + 1 <= j else 0) for k in range(R[i, j - 1], R[i + 1, j] + 1)]

            minWeight = min(weights)
            minIndex = weights.index(minWeight) + R[i, j - 1]

            R[i, j] = minIndex
            P[i, j] = minWeight + W[i, j]

    return R


def parentsFromRoots(R):
    """
    Build the parents array from the roots matrix
    For each node, it gives the index of its parent in the optimal
    subtree, with the root of the subtree being its own parent

    :param R: the matrix of the roots of the optimal subtrees
    """
    parents = [-1] * len(R)
    buildSubTree(R, parents, -1, 0, len(R) - 1)

    # Locate the root and replace it with its own index
    generalRoot = parents.index(-1)
    parents[generalRoot] = generalRoot

    return parents


def buildSubTree(R, parents, calledFrom, i, j):
    """
    Recursively build the array of parents

    :param R: the matrix of the root of the optimal subtrees 
    :param parents: the array to be built, containing the list of the parents of each node
    :param calledFrom: the index of the node from which this function was called
    :param i: index of the left bound of the subtree
    :param j: index of the right bound of the subtree
    """
    parents[R[i, j]] = calledFrom
    if not (i <= R[i, j] and R[i, j] <= j):
        raise BaseException("The root is out of the bounds of the subtree")
    if i < j:
        if i <= R[i, j] - 1:
            # In case we indeed have a left subtree
            buildSubTree(R, parents, R[i, j], i, R[i, j] - 1)
        if R[i, j] + 1 <= j:
            # In case we indeed have a right subtree
            buildSubTree(R, parents, R[i, j], R[i, j] + 1, j)
