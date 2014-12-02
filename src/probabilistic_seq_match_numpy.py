# -*- coding: utf-8 -*-
#
# Code containing functions for calculatinf sequence match probabilities
#
import timeit
import numpy as np
import random
import sys

def main():
    base_prob_precompute = {}
    for score in range(1, 51):
        base_prob_precompute[score] = base_prob(score)
    
    seqa = "ANTCATCTTCAACTCGGCAGGAACCCACAAGTGTGCATGTGTGGCTCGGAGGCTTCAGCTGGGGCCCCGN"
    quala = [random.randint(1, 50) for i in range(len(seqa))]
    seqb = seqa[::-1]
    qualb = quala[::-1]

    # Wrap the non numpy and test
    num = 1
    print("Pure python: {0} times".format(num))
    py_wrapped = wrapper(sequences_match_prob, seqa, quala, seqb, qualb,
                         base_prob_precompute, 0)
    res = timeit.timeit(py_wrapped, number=num)
    print(res)

    # Make arrays
    seqa_np = np.array(list(seqa))
    seqb_np = np.array(list(seqb))
    quala_np = np.array(list(quala))
    qualb_np = np.array(list(qualb))

    # Test numpy
    print("Numpy: {0} times".format(num))
    np_wrapped = wrapper(numpy_seq_prob, seqa_np, quala_np, seqb_np, qualb_np,
                         base_prob_precompute)
    res = timeit.timeit(np_wrapped, number=num)
    print(res)

def numpy_seq_prob(seqa, quala, seqb, qualb, base_prob_precompute):

    # Convert quals to base probabilites
    proba = np.power(10, -quala/10)
    probb = np.power(10, -qualb/10)

    # Result array
    match_prob = np.empty(proba.shape)
    match_prob.fill(np.nan)

    # Where either base is N, set match prob to 0.25
    n = np.logical_or(seqa == "N", seqb == "N")
    match_prob[n] = 0.25
    
    # If bases match calc prob
    match = np.logical_and(seqa == seqb, np.logical_not(n))
    match_prob[match] = ((1 - proba[match]) * (1 - probb[match]) +
                         (proba[match] * probb[match]) / 3)

    # If bases don't match, calc prob
    no_match = np.logical_and(np.logical_not(n), np.logical_not(match))
    match_prob[no_match] = ((1 - proba[no_match]) * probb[no_match] / 3 +
                            (1 - probb[no_match]) * proba[no_match] / 3 +
                            2 * proba[no_match] * probb[no_match] / 9)

    print(np.prod(match_prob))
    print(match_prob)
    return np.prod(match_prob)

def wrapper(func, *args, **kwargs):
    def wrapped():
        return func(*args, **kwargs)
    return wrapped

def sequences_match_prob(a_seq, a_qual, b_seq, b_qual, base_prob_precompute,
                         stop_thresh):
    """ Given two sequences and their quality scores
    """

    # Calc prob
    prob = 1.0

    # debug
    match_probs = []

    # For each base
    for i in range(len(a_seq)):
        # If Either is N, then prob is 0.25
        if a_seq[i] == 'N' or b_seq[i] == 'N':
            match_prob = 0.25
        else:
            # Calculate the base probabilities
            #prob_A = base_prob_precompute[a_qual[i]]
            #prob_B = base_prob_precompute[b_qual[i]]
            # Calculate the base probabilities
            prob_A = base_prob(a_qual[i])
            prob_B = base_prob(b_qual[i])
            # Calc probability of a match
            if a_seq[i] == b_seq[i]:
                match_prob = match_given_match_prob(prob_A, prob_B)
            elif a_seq[i] != b_seq[i]:
                match_prob = match_given_mismatch_prob(prob_A, prob_B)
        # Combine with overall prob
        prob = prob * match_prob
        # DEBUG
        match_probs.append(match_prob)
        # Return if less than threshold
        if prob < stop_thresh:
            return prob
    print(prob)
    print(match_probs)
    return prob

def match_given_mismatch_prob(x_prob, y_prob):
    """ Gives the prob of true match given two bases match.
    """
    return (  (1.0/3) * (1 - x_prob) * y_prob
            + (1.0/3) * (1 - y_prob) * x_prob
            + (2.0/9) * x_prob * y_prob )

def match_given_match_prob(x_prob, y_prob):
    """ Gives the prob of true match given two bases match.
    """
    return (1 - x_prob) * (1 - y_prob) + (x_prob * y_prob) / 3

def base_prob(phred_score):
    """ Returns the probabilty that a base is incorrect, given its
        Phred score.
    """
    prob = 10.0**(-float(phred_score)/10)
    return prob

def phred_score(letter, offset, ascii):
    """ Returns the Phred score of a fastq quality character.
    """
    # Offset the string
    offset_string = ascii[offset - 33:]
    # Find the characters score
    score = 0
    while score < len(offset_string):
        if letter == offset_string[score]:
            return score
        score += 1
    # If no score is found then there must be an error
    raise ValueError

if __name__ == '__main__':
    main()
