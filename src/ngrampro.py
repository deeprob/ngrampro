import numpy as np


def check_sequence(sequence):
    """
    A function to verify that there are no unknown amino acids in the protein sequence provided.
    :param sequence: The protein sequence.
    :return:
    """
    AA = 'ACDEFGHIKLMNPQRSTVWY'
    aa_set = set(sequence)
    for aa in aa_set:
        if aa not in AA:
            raise ValueError("Unknown Amino Acid detected in protein sequence.")
    return


def check_sequences(sequences):
    """
    A function to check the given list of protein sequences comply with two rules:
    1. The string that represents the protein sequence is uppercase
    2. The protein sequence contains only valid amino acids.
    :return:
    """
    if all([seq.isupper() for seq in sequences]):
        list(map(check_sequence, sequences))
    else:
        raise ValueError("Protein sequence has both upper and lower case strings. Please convert it to upper case "
                         "characters only.")
    return


class GAA:
    """
    A class that changes a given list of amino acid sequences to their grouped representation. Grouping is achieved by
    categorizing the amino acids into five different classes based on their physicochemical properties. The Amino Acids
    are grouped into five classes: aliphatic group (G,A,V,L,M,I), aromatic group (F,Y,W), positive charge group (K,R,H),
    negatively charged group (D,E) and uncharged group (S,T,C,P,N,Q).

    Attributes
    ----------
    groups : dictionary
        A dictionary that maps an amino acid to its group.

    Methods
    -------
    transform_map(sequence):
        Converts a single protein sequence to its grouped representation.

    transform(sequences):
        Converts multiple protein sequences to their grouped representation.
    """
    def __init__(self):
        """
        Constructs necessary attributes for the GAA object.
        """
        self.groups = {'G': 'a', 'A': 'a', 'V': 'a', 'L': 'a', 'M': 'a', 'I': 'a', 'F': 'b', 'Y': 'b', 'W': 'b',
                       'K': 'c', 'R': 'c', 'H': 'c', 'D': 'd', 'E': 'd', 'S': 'e', 'T': 'e', 'C': 'e', 'P': 'e',
                       'N': 'e', 'Q': 'e'}

    def transform_map(self, sequence):
        """
        Converts a single protein sequence to its grouped representation.
        :param sequence: A string that represents the protein sequence.
        :return: A string that represents the grouped representation of the protein sequence.
        """
        check_sequence(sequence)
        return ''.join([self.groups[aa] for aa in sequence])

    def transform(self, sequences):
        """Transforms an a list of sequences to a grouped representation"""
        X_gaa = list(map(self.transform_map, sequences))
        return np.array(X_gaa)


class Ngram:
    """
    A class used to represent a sequence as a one-hot encoded representation of its n-grams.

    Attributes
    ----------
    sequences : list like representation of str
        list of protein sequences.
    n : int
        the k-mer size.
    step : int
        the step size to consider while iterating over the sequences.
    inc_count : bool
        A boolean argument whether to keep the counts of the n-grams present in the sequence.
    encoder_dict : dict
        A dictionary that stores the identified n-grams in the given list of sequences as keys and their counts as
        values.
    keep_v : int
        Only n-grams that are present beyond this number on the given list of sequences will be stored in the
        encoder dict.
    Methods
    -------
    fit():
        Creates the encoder dict from the list of sequences.

    transform(sequences):
        Given a list of sequences, this function creates a one-hot encoded n-gram representation of those sequences.
        It is necessary to use fit first before using transform since the encoder dict needs to be created.
    """
    def __init__(self, sequences, n, step=1, inc_count=False, keep_v=1):
        """
        Constructs necessary attributes for the Ngram object.

        Parameters
        ----------
        sequences : list like representation of str
            list of protein or grouped protein sequences or any sequences in general.
        n : int
            the n-gram size.
        step : int
            the step size to consider while iterating over the sequences.
        inc_count : bool
            A boolean argument whether to keep the counts of the k-mers present in the sequence.
        keep_v : int
            Only n-grams that are present beyond this number on the given list of sequences will be stored in the
            encoder dict.
        """
        self.sequences = sequences
        self.n = n
        self.step = step
        self.inc_count = inc_count
        self.keep_v = keep_v
        self.encoder_dict = {}

    def fit(self):
        """
        Creates the encoder dict from the list of sequences. Encoder dict stores the identified n-grams as keys and
        their counts (if inc_count is True else 1) as values
        :return: None
        """
        ngram_dict = dict()
        for sequence in self.sequences:
            i = 0
            while i + self.n <= len(sequence):
                if type(sequence) == str:
                    seq = sequence[i:i + self.n]
                else:
                    seq = ''.join(sequence[i:i + self.n])

                if seq in ngram_dict:
                    ngram_dict[seq] += 1
                else:
                    ngram_dict[seq] = 1
                i += self.step

        feature_list = sorted([k for k, v in ngram_dict.items() if v >= self.keep_v])
        self.encoder_dict = dict(zip(feature_list, list(range(len(feature_list)))))
        return

    def transform(self, sequences):
        """
        Given a list of sequences, this function creates a one-hot encoded n-gram representation of those sequences. It
        is necessary to use fit first before using transform since the encoder dict needs to be created.
        :param sequences: a list of sequences which needs to be transformed into a one-hot encoded representation of
        n-grams.
        :return: an array of one-hot encoded protein sequences or any other sequences.
        """
        if not self.encoder_dict:
            raise ValueError('Need to fit first')

        X_motif = []

        for sequence in sequences:
            ind_vector = np.zeros(len(self.encoder_dict))

            i = 0
            while i + self.n <= len(sequence):
                if type(sequence) == str:
                    seq = sequence[i:i + self.n]
                else:
                    seq = ''.join(sequence[i:i + self.n])

                if seq in self.encoder_dict:
                    if self.inc_count:
                        ind_vector[self.encoder_dict[seq]] += 1
                    else:
                        ind_vector[self.encoder_dict[seq]] = 1

                i += self.step

            X_motif.append(ind_vector)

        return np.array(X_motif)


class NGModel:
    """
    A class that creates a one-hot encoded representation of protein sequences from a given training, validation and
    optional test set. The one-hot encoded representation is created from k-mers present in the training set. It is
    ensured that there is no data leakage from training to validation or test set.

    Attributes
    ----------
    x_train_raw : list like representation of str
        list of protein sequences from the given training set.
    x_valid_raw : list like representation of str
        list of protein sequences from the given validation set.
    ng : object
        the Ngram object.
    x_train : array
        one-hot encoded representation of protein sequences in the training set.
    x_valid : array
        one-hot encoded representation of protein sequences in the validation set.
    x_test : array or None
        one-hot encoded representation of protein sequences in the optional test set.
    """
    def __init__(self, x_train, x_valid, x_test=None, k=3, s=1, inc_count=False, keep_v=1):
        """
        Constructs the necessary attributes for the NGModel object

        Parameters
        ----------
        x_train : list like representation of str
            list of protein sequences from the given training set.
        x_valid : list like representation of str
            list of protein sequences from the given validation set.
        x_test : list like representation of str (optional)
            list of protein sequences from the given validation set.
        k : int
            the k-mer size.
        s : int
            the step size to consider while iterating over the sequences.
        inc_count : bool
            A boolean argument whether to keep the counts of the k-mers present in the sequence.
        keep_v : int
            Only k-mers that are present beyond this number on the given list of sequences will be stored in the
            encoder dict.
        """
        check_sequences(x_train), check_sequences(x_valid)
        self.x_train_raw, self.x_valid_raw = x_train, x_valid
        self.ng = Ngram(self.x_train_raw, k, s, inc_count, keep_v)
        self.ng.fit()
        self.x_train, self.x_valid = self.ng.transform(self.x_train_raw), self.ng.transform(self.x_valid_raw)

        if x_test is not None:
            check_sequences(x_test)
            self.x_test = self.ng.transform(x_test)
        else:
            self.x_test = None


class GAANGModel:
    """
    A class that creates a one-hot encoded representation of grouped protein sequences from a given training, validation
    and optional test set. The one-hot encoded representation is created from k-mers present in the training set by
    first transforming the sequences into a grouped representation and then applying Ngram. It is
    ensured that there is no data leakage from training to validation or test set.

    Attributes
    ----------
    x_gaa_train : array of str
        array of grouped protein sequences from the given training set.
    x_gaa_valid : array of str
        array of grouped protein sequences from the given validation set.
    x_gaa_test : array of str or None
        array of grouped protein sequences from the given validation set.
    ga : object
        the GAA object.
    ng : object
        the Ngram object.
    x_train : array
        one-hot encoded representation of protein sequences in the training set.
    x_valid : array
        one-hot encoded representation of protein sequences in the validation set.
    x_test : array or None
        one-hot encoded representation of protein sequences in the optional test set.
    """
    def __init__(self, x_train, x_valid, x_test=None, k=3, s=1, inc_count=False, keep_v=1):
        """
        Constructs the necessary attributes for the GAANGModel object

        Parameters
        ----------
        x_train : list like representation of str
            list of protein sequences from the given training set.
        x_valid : list like representation of str
            list of protein sequences from the given validation set.
        x_test : list like representation of str (optional)
            list of protein sequences from the given validation set.
        k : int
            the k-mer size.
        s : int
            the step size to consider while iterating over the sequences.
        inc_count : bool
            A boolean argument whether to keep the counts of the k-mers present in the sequence.
        keep_v : int
            Only k-mers that are present beyond this number on the given list of sequences will be stored in the
            encoder dict.
        """

        check_sequences(x_train), check_sequences(x_valid)
        self.ga = GAA()
        self.x_gaa_train = self.ga.transform(x_train)
        self.x_gaa_valid = self.ga.transform(x_valid)
        self.ng = Ngram(self.x_gaa_train, k, s, inc_count, keep_v)
        self.ng.fit()
        self.x_train, self.x_valid = self.ng.transform(self.x_gaa_train), self.ng.transform(self.x_gaa_valid)

        if x_test is not None:
            check_sequences(x_test)
            self.x_gaa_test = self.ga.transform(x_test)
            self.x_test = self.ng.transform(self.x_gaa_test)
        else:
            self.x_gaa_test = None
            self.x_test = None
