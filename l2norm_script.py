import statistics
import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

import minhash
import sequences
import align

mh = minhash.SequenceComparison()
random.seed(0)

# We want to compare the Levenshtein similarity metric with
# the euclidean distance between k-mer frequency vectors for
# of strings. The k-mer count representation implementation
# can be found in the minhash module.

# We will compare these two similarity metrics

seq_len = 1000
s1 = s2 = sequences.random_sequence(seq_len)
num_edits = int(seq_len * 0.4)

k_range = [2, 4, 6, 8, 10]
data = pd.DataFrame()

for e in range(num_edits):
    while mh.edit_dist(s1, s2) <= e:
        s2 = sequences.random_edits(s2, 1)
    for k in k_range:
        est_dist = mh.kmer_vec_dist(s1, s2, k)
        data.at[e, f'k={k}'] = est_dist

print(data.head())

def plot():
    fontdict = {'fontsize': 15}
    plt.figure(figsize=(12, 6))
    for column in data.columns:
        sns.regplot(x=data.index, y=column, label=column,
            data=data,
            scatter_kws={'s': 50, 'label': column})
    plt.xlabel('Levenshtein Distance', fontdict=fontdict)
    plt.ylabel('Euclidean Distance Between k-mer Frequency Vectors', fontdict=fontdict)
    plt.title(f'Figure 1: Euclidean Distance. String length: {seq_len}.', fontdict=fontdict)
    plt.legend()
    plt.grid(True)
    plt.show()

plot()
