# PolyaModelsYeast

This repository contains the necessary scripts to make predictions using the *S. cerevisiae*-specific models **PolyaClassifier**, **PolyaCleavage**, and **PolyaStrength**, convolutional neural network models that predict the classification, cleavage site distribution, and strength of a polyA site, respectively. We simultaneously developed a *S. pombe*-specific **PolyaClassifier** model. We developed this suite of models to delineate the *cis*-regulatory elements determining polyA site selection and usage in yeast.

Contact *zhe.ji (at) northwestern.edu* with any questions.


## Making New Predictions

### Running PolyaModels for yeast requires the following packages be installed

- Python == 3.6
- Tensorflow == 2.1.0
- Keras == 2.3.1
- NumPy == 1.19.1
- Pandas == 1.1.5
- pyfaidx == 0.5.9
- Isolearn

### Usage

PolyaClassifier, PolyaCleavage, and PolyaStrength can be used to make new predictions from a genomic location or sequence. Genomic locations must be given as a string with the format "chrom:position:strand" and the reference genome FASTA and chromosome sizes files must be provided.

**predictor_tool/yeast_prediction.py**
> This file contains the predictor tool to make new predictions. It is designed to be used as a command-line tool, which users can invoke as shown in the examples below.

**predictor_tool/yeast_utilities.py**
> This file contains the helper functions used by the predictor tool. For example, functions to extract the genomic sequence if needed, build and reload the PolyaModels, and flow batches of data into the models for predictions are here.

**predictor_tool/S_cerevisiae.PolyaClassifier.h5**
> The trained model weights for the *S. cerevisiae*-specific PolyaClassifier model.

**predictor_tool/S_cerevisiae.PolyaCleavage.h5**
> The trained model weights for the *S. cerevisiae*-specific PolyaCleavage model.

**predictor_tool/S_cerevisiae.PolyaStrength.h5**
> The trained model weights for the *S. cerevisiae*-specific PolyaStrength model.

**predictor_tool/S_pombe.PolyaClassifier.h5**
> The trained model weights for the *S. pombe*-specific PolyaClassifier model.


### Example prediction from sequence

```sh
python yeast_prediction.py from_sequence -m 'S.cer PolyaClassifier' -s 'CGTGTGTTTTTGTTCTCTTATAACCGAGCTGCTTACTTATTATTATTTCACCTTCTCTTTTTATTTATACTTATAATTATTTATTCTTTACATACTGTTACAAGAAACTCTTTTCTACATTAATTGCATAAAGTGTCAATCAGCACATCCTCTACATCGCTATCAACAACAAATTTGACAAACCTGCCTATATCTTCAGGAACGACTGCTGCATCGCTACCACCACTACTTGTGAAGTCC'
```

This will give the following output: 

```
Sequence: CGTGTGTTTTTGTTCTCTTATAACCGAGCTGCTTACTTATTATTATTTCACCTTCTCTTTTTATTTATACTTATAATTATTTATTCTTTACATACTGTTACAAGAAACTCTTTTCTACATTAATTGCATAAAGTGTCAATCAGCACATCCTCTACATCGCTATCAACAACAAATTTGACAAACCTGCCTATATCTTCAGGAACGACTGCTGCATCGCTACCACCACTACTTGTGAAGTCC
PolyaClassifier probability: 0.999801
```

A genomic location with this sequence is expected to be a polyA cleavage site with a classification probability >0.999 from PolyaClassifier.

### Example prediction from genomic location

**Note:** Predictions can be made from genomic locations but the genome FASTA and chrom.sizes files need to be provided by the user.

```sh
python yeast_prediction.py from_position -m 'S.cer PolyaClassifier' -p 'I:42717:-' -g './genome.fa' -c './chrom.sizes' 
```

This will give the following output:

```
Sequence: CGTGTGTTTTTGTTCTCTTATAACCGAGCTGCTTACTTATTATTATTTCACCTTCTCTTTTTATTTATACTTATAATTATTTATTCTTTACATACTGTTACAAGAAACTCTTTTCTACATTAATTGCATAAAGTGTCAATCAGCACATCCTCTACATCGCTATCAACAACAAATTTGACAAACCTGCCTATATCTTCAGGAACGACTGCTGCATCGCTACCACCACTACTTGTGAAGTCC
PolyaClassifier probability: 0.999801
```

At chromosome I, position 42,717 on the reverse strand, a putative polyA site at this location is expected with a classification probability >0.999 from PolyaID.

