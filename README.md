# ngrampro
A python tool to generate k-mers from protein sequences based on N-grams


# ngrampro module main classes
There are two main classes that can be used to generate one-hot encoded k-mer representation of protein sequences. 
They are:
1. NGModel
2. GAANGModel

## NGModel
NGModel can be used when the user wants to create a one-hot encoded k-mer representation of protein sequences. It 
requires training and validation set of protein sequences. An optional test set may also be provided. It first creates 
an encoder dict that maps all the k-mers present in the training set of sequences to the number of times they were found
 in the training set. The one-hot encoded representation for the training, validation and test set is created from this 
 encoder dict. Please note that neither the validation nor the test set is used to create the encoder dict thus ensuring
no data leakage from training to validation or test.

## GAANGModel
GAANGModel can be used when the user wants to create a one-hot encoded k-mer representation of grouped protein 
sequences. It requires training and validation set of protein sequences. An optional test set may also be provided. 
It is similar in operation NGModel with an additional preprocessing step. In the preprocessing step, the amino acids in 
the protein sequences are first categorized into one of the five pre-defined groups based on their physicochemical 
properties. Then for each protein sequence, a new grouped protein representation is created where an amino acid is 
represented by its group. The grouped protein representation is used to create the one-hot encoded k-mer based protein 
sequence representation.


# ngrampro module sub classes
There are two subclasses in the module which are used by the main classes.
1. GAA
2. Ngram

## GAA
GAA class can be used to create a grouped protein sequence representation from the original sequence.

## Ngram
Ngram class can be used to create a one-hot encoded representation of sequences based on the n-grams present in those 
sequences. This is the only class in this module that can work with any sequences, not just protein sequences.
