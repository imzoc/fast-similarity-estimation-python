import random

class Sequences:
    def __init__(self) -> None:
        self.alphabet = ['A', 'C', 'T', 'G']

    def random_sequence(self, l):
        return ''.join(random.choices(self.alphabet, k=l))
    def random_sequences(self, n, l):
        return tuple(self.random_sequence(l) for _ in range(n))
        
def random_sequences(n, l):
    s = Sequences()
    return s.random_sequences(n, l)

def random_sequence(l):
    s = Sequences()
    return s.random_sequence(l)

def random_edits(seq, n):
    seq = list(seq)
    alphabet = ['A', 'C', 'T', 'G']
    for _ in range(n):
        index = random.randint(0, len(seq) - 1)
        kind = random.choice(['substitution', 'insertion', 'deletion'])
        if kind == 'substitution':
            seq[index] = random.choice(alphabet)
        if kind == 'insertion':
            seq.insert(index, random.choice(alphabet))
        if kind == 'deletion':
            seq.pop(index)
    return ''.join(seq)

