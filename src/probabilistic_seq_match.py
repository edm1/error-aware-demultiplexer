# -*- coding: utf-8 -*-
#
# Code containing functions for calculatinf sequence match probabilities
#

def sequences_match_prob(a_seq, a_qual, b_seq, b_qual, base_prob_precompute,
                         stop_thresh):
    """ Given two sequences and their quality scores
    """

    # Calc prob
    prob = 1.0

    # For each base
    for i in range(len(a_seq)):
        # If Either is N, then prob is 0.25
        if a_seq[i] == 'N' or b_seq == 'N':
            match_prob = 0.25
        else:
            # Calculate the base probabilities
            prob_A = base_prob_precompute[a_qual[i]]
            prob_B = base_prob_precompute[b_qual[i]]
            # Calc probability of a match
            if a_seq[i] == b_seq[i]:
                match_prob = match_given_match_prob(prob_A, prob_B)
            elif a_seq[i] != b_seq[i]:
                match_prob = match_given_mismatch_prob(prob_A, prob_B)
        # Combine with overall prob
        prob = prob * match_prob
        # Return if less than threshold
        if prob < stop_thresh:
            return prob
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
