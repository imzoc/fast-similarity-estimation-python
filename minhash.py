import math
import editdistance
import random

class SequenceComparison:
    """ This class contains methods relating to string
    similarity estimation methods.
    """
    def __init__(self) -> None:
        self.alphabet = ['A', 'C', 'T', 'G']

    def kmer_space(self, seq, k=5):
        ''' This method creates and returns a tensor
        representing the k-mer content of the sequence.
        It uses class KmerSpace.
        '''
        T = KmerSpace(k)
        
        # sliding window
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            indices = [self.alphabet.index(c) for c in kmer]
            T.update(indices)

        return T
    
    def random_hash_function(self, seed):
        ''' this function returns reproducibly random hash
        functions using a seed
        '''
        return lambda s: abs(hash((s, seed))) % 10000 + 1
    
    def edit_dist(self, s1, s2):
        ''' returns the edit distance between two strings
        (implemented in a known library)
        '''
        return editdistance.eval(s1, s2)

    def kmer_jaccard_sim(self, s1, s2, k):
        ''' returns the euclidean distance between the kmer
        vector representations of strings.
        '''
        s1_vec = self.kmer_space(s1, k=k).as_vector()
        s2_vec = self.kmer_space(s2, k=k).as_vector()
        return self.jaccard_sim(s1_vec, s2_vec)

    def kmer_vec_dist(self, s1, s2, k):
        ''' returns the euclidean distance between the kmer
        vector representations of strings.
        '''
        s1_vec = self.kmer_space(s1, k=k).as_vector()
        s2_vec = self.kmer_space(s2, k=k).as_vector()
        return self.l2norm(s1_vec, s2_vec)
    
    def kmer_vec_cosine_sim(self, s1, s2, k):
        ''' returns the euclidean distance between the kmer
        vector representations of strings.
        '''
        s1_vec = self.kmer_space(s1, k=k).as_vector()
        s2_vec = self.kmer_space(s2, k=k).as_vector()
        return self.cosine_similarity(s1_vec, s2_vec)

    def minimizer_comparison(self, s1, s2, w, k, f):
        ''' sketches each string using minimizers, then
        compares each minimizer. Returns the fraction of
        minimizers that are the same.
        '''
        s1_minimizers = []
        s2_minimizers = []
        for i in range(len(s1) - w + 1):
            substring = s1[i:i+w]
            kmers = self.kmer_set(substring, k)
            minimizer = min(kmers, key=f)
            s1_minimizers.append(minimizer)

        for i in range(len(s2) - w + 1):
            substring = s2[i:i+w]
            kmers = self.kmer_set(substring, k)
            minimizer = min(kmers, key=f)
            s2_minimizers.append(minimizer)

        return sum(x1==x2 for x1, x2 in zip(s1_minimizers, s2_minimizers))

    def kmer_set(self, s, k):
        ''' return the set of kmers over s '''
        return [s[i:i+k] for i in range(len(s) + k - 1)]

    def l2norm(self, v1, v2):
        ''' returns the l2 norm of two vectors
        (the euclidean distance between them)
        '''
        return math.sqrt(
            sum([(x1-x2)**2 for x1, x2 in zip(v1, v2)])
        )

    def cosine_similarity(self, v1, v2):
        ''' returns the cosine similarity between two vectors
        '''
        dot_product = sum(x1 * x2 for x1, x2 in zip(v1, v2))
        magnitude_v1 = math.sqrt(sum(x ** 2 for x in v1))
        magnitude_v2 = math.sqrt(sum(x ** 2 for x in v2))
        
        if magnitude_v1 == 0 or magnitude_v2 == 0:
            return 0.0
        return dot_product / (magnitude_v1 * magnitude_v2)

    def jaccard_sim(self, v1, v2):
        ''' returns the jaccard similarity score between
        two lists of kmer frequencies
        '''
        intersection = sum([min(x1, x2) for x1, x2 in zip(v1, v2)])
        union = sum([max(x1, x2) for x1, x2, in zip(v1, v2)])
        return intersection / union

class KmerSpace:
    ''' This class contains a flat tensor space. Each k-mer is a 
    vector in this space, with each index representing a 
    dimension. A full string is represented by the entire
    tensor space.
    '''
    def __init__(self, k) -> None:
        self.k = k
        self.T = [0] * (4 ** k)
        self.total = 0

    def __str__(self):
        return str(self.T)
    
    def __getitem__(self, indices):
        index = self._get_flat_index(indices)
        return self.T[index]

    def _get_flat_index(self, indices):
        index = 0
        for i, val in enumerate(indices):
            index += val * (4 ** (self.k - i - 1))
        return index

    def __setitem__(self, indices, item):
        index = self._get_flat_index(indices)
        self.T[index] = item

    def update(self, indices):
        self.total += 1
        index = self._get_flat_index(indices)
        self.T[index] += 1

    def as_vector(self):
        return self.T



