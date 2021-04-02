# ngrampro (ngram for protein sequences)
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

```python
!pip install ngrampro
```

    Collecting ngrampro
      Using cached ngrampro-0.0.5-py3-none-any.whl (5.3 kB)
    Requirement already satisfied: numpy in /Users/{user_name}/opt/anaconda3/lib/python3.7/site-packages (from ngrampro) (1.17.2)
    Installing collected packages: ngrampro
    Successfully installed ngrampro-0.0.5

# ngrampro usage

```python
import ngrampro as npro
```

## Defining some random protein sequences
This will be used for illustrating how ngrampro can be used


```python
import random
random.seed(0)

AA = 'ACDEFGHIKLMNPQRSTVWY'
train_sequences = ["".join([random.choice(AA) for _ in range(3)]) for _ in range(5)]
valid_sequences = ["".join([random.choice(AA) for _ in range(3)]) for _ in range(5)]
test_sequences = ["".join([random.choice(AA) for _ in range(3)]) for _ in range(5)]
print(train_sequences, valid_sequences, test_sequences)
```

    ['PQC', 'KTS', 'PLS', 'NWH', 'TFL'] ['FEY', 'KVY', 'FLE', 'DMS', 'VEN'] ['QMY', 'HVS', 'RTK', 'CVA', 'DPA']


## Using the NGModel object

NGModel object creates an one-hot encoded representation of the k-mers present in protein sequences using the ngram modeling principle. It requires two sets of protein sequences termed training and validation and an optional test set of protein sequences. The protein sequences can be a entered as a list like object of strings. The user can enter the k-mer size they would like as well as the minimum k-mer count among the k-mers found in the training set that would be used to create the one-hot encoded representation using the k and keep_v arguments respectively. Thus if the minimum k-mer count is 1, only k-mers found at least once in the training set will be used to create the one-hot encoded dict in the training, validation and test set.  


```python
# The ngmodel class with create an encoder dict that maps all the k-mers found in the training set of protein 
# sequences to their corresponding count. This encoder dict is used to subsequently create a one-hot encoded
# representation of proteins sequences in the training, validation and test set. 

ngmodel = npro.NGModel(train_sequences, valid_sequences, test_sequences, k=2, keep_v=0)
```


```python
# The attributes x_train, x_valid and x_test can be used to recover the one hot encoded representation of 
# protein sequences in the training, validation and test set developed by ngmodel 
ngmodel.x_train
```




    array([[0., 0., 0., 0., 0., 1., 1., 0., 0., 0.],
           [0., 1., 0., 0., 0., 0., 0., 0., 1., 0.],
           [0., 0., 1., 0., 1., 0., 0., 0., 0., 0.],
           [0., 0., 0., 1., 0., 0., 0., 0., 0., 1.],
           [1., 0., 0., 0., 0., 0., 0., 1., 0., 0.]])




```python
ngmodel.x_valid
```




    array([[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
           [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
           [1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
           [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
           [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]])




```python
ngmodel.x_test
```




    array([[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
           [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
           [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
           [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
           [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]])



## Using the GAANGModel object 

GAANGModel works exactly like NGModel with the same required and optional arguments. The only difference is in the way it creates the one-hot encoded representation under the hood. It has an additional preprocessing step that creates a grouped representation of the protein sequences where each amino acid is replaced by its designated group/category based on its physicochemical property. GAANGModel can be more useful when the number of protein sequences in the training set is low.


```python
gaangmodel = npro.GAANGModel(train_sequences, valid_sequences, test_sequences, k=2, keep_v=0)
```


```python
gaangmodel.x_train
```




    array([[0., 0., 0., 0., 0., 0., 1.],
           [0., 0., 0., 1., 0., 0., 1.],
           [1., 0., 0., 0., 1., 0., 0.],
           [0., 0., 1., 0., 0., 1., 0.],
           [0., 1., 0., 0., 0., 1., 0.]])




```python
gaangmodel.x_valid
```




    array([[0., 0., 0., 0., 0., 0., 0.],
           [0., 0., 0., 0., 0., 0., 0.],
           [0., 1., 0., 0., 0., 0., 0.],
           [1., 0., 0., 0., 0., 0., 0.],
           [0., 0., 0., 0., 0., 0., 0.]])




```python
gaangmodel.x_test
```




    array([[0., 0., 0., 0., 1., 0., 0.],
           [1., 0., 0., 0., 0., 0., 0.],
           [0., 0., 0., 1., 0., 0., 0.],
           [0., 0., 0., 0., 1., 0., 0.],
           [0., 0., 0., 0., 1., 0., 0.]])
