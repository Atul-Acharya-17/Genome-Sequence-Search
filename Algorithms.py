import time
import pickle

"""
    This is the python file that contains all the important algorithms for our project
    The brute_force_search, knuth_morris_pratt and rabin_karp search functions can be found below
    
    All of these functions keep track of the time elapsed and the number of times the query appears in the genome
    sequence. It also keeps track of the number of character comparisons made by these 3 algorithms and stores 
    the length of the genome sequence, length of the query and the number of character comparisons in a binary file 
"""


# Implementation of the brute force algorithm
def brute_force_search(sequence, pattern):
    comparisons = 0
    positions = []
    count = 0
    start = time.time()
    for i in range(len(sequence) - len(pattern) + 1):
        t = 0
        for j in range(len(pattern)):
            comparisons += 1
            if sequence[i + j] != pattern[j]:
                break
            else:

                t = t + 1
        if t == len(pattern):
            count += 1
            positions.append(i + 1)

    end = time.time()
    time_elapsed = end - start

    file = "brute_force_observations.pkl"
    dump_observations(file, sequence, pattern, comparisons)

    return positions, count, time_elapsed, comparisons


# Implementation of the Knuth Morris Pratt algorithm
def knuth_morris_pratt(sequence, pattern):
    comparisons = 0
    positions = []
    i = 0
    j = 0
    count = 0
    start = time.time()
    pi_table = kmp_preprocess(pattern)
    while i < len(sequence):
        comparisons += 1
        if sequence[i] == pattern[j]:
            i += 1
            j += 1
        if j == len(pattern):
            positions.append((i - j + 1))
            j = pi_table[j - 1]
            count += 1
        elif i < len(sequence) and sequence[i] != pattern[j]:
            if j == 0:
                i += 1
            else:
                j = pi_table[j - 1]

    end = time.time()
    time_elapsed = end - start

    file = "knuth_morris_pratt_observations.pkl"
    dump_observations(file, sequence, pattern, comparisons)

    return positions, count, time_elapsed, comparisons


# This function calculates the pi table which is used in the KMP Algorithm
def kmp_preprocess(pattern):
    i = 1
    j = 0
    pi_table = [0] * len(pattern)

    while i < len(pattern):
        if pattern[i] == pattern[j]:
            j += 1
            pi_table[i] = j
            i += 1
        else:
            if j == 0:
                pi_table[i] = 0
                i += 1
            else:
                j = pi_table[j - 1]

    return pi_table


# Implementation of the Rabin Karp algorithm
def rabin_karp(sequence, pattern):
    comparisons = 0
    n = len(sequence)
    m = len(pattern)
    q = 101
    base = 5
    x = 1
    index = []
    hash_q = 0
    hash_g = 0
    start = time.time()
    for i in range(0, m - 1):
        x = (x * base) % q
    for i in range(0, m):
        hash_g = (base * hash_g + ord(sequence[i])) % q
        hash_q = (base * hash_q + ord(pattern[i])) % q
    for i in range(0, (n - m + 1)):
        comparisons += 1
        if hash_q == hash_g:
            for j in range(0, m):
                comparisons += 1
                if pattern[j] != sequence[i + j]:
                    break

            if j == m - 1:
                index.append(i + 1)
        if i < n - m:
            hash_g = (base * (hash_g - ord(sequence[i]) * x) + ord(sequence[i + m])) % q
            if hash_g < 0:
                hash_g = hash_g + q

    end = time.time()
    time_elapsed = end - start

    file = "rabin_karp_observations.pkl"
    dump_observations(file, sequence, pattern, comparisons)

    return index, len(index), time_elapsed, comparisons


# This function loads all our previous observations into the program
def load_observations(filename):
    with open(filename, "rb") as f:
        observations = pickle.load(f)
    f.close()
    return observations


# This function appends the records with the new record
def dump_observations(file, sequence, pattern, comparisons):
    try:
        observations = load_observations(file)
        observations.append((len(sequence), len(pattern), comparisons))
    except Exception as e:
        print(e)
        observations = [(len(sequence), len(pattern), comparisons)]

    with open(file, 'wb') as f:
        pickle.dump(list(observations), f)

    f.close()
