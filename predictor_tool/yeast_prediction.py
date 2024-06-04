#!/usr/bin/env python
# coding: utf-8


##############################
## IMPORTS
##############################


import yeast_utilities as utils

import os, argparse
import numpy as np

from contextlib import redirect_stderr

with redirect_stderr(open(os.devnull, "w")):
	from keras.models import load_model



##############################
## HELPER FUNCTIONS
##############################

class Predictor(object):

	def __init__(self, modelpath):
		self.model = load_model(modelpath)

	def predict(self, encoded_input):
		return self.model.predict(encoded_input)



class Ensemble(object):

	def __init__(self, modelpaths):
		self.models = [Predictor(mp) for mp in modelpaths]

	def ind_predict(self, encoded_input):
		return [m.predict(encoded_input) for m in self.models]

	def predict(self, encoded_input, agg_function = np.mean):
		return agg_function(self.ind_predict(encoded_input), axis = 0)



def report_predictions(lmodel, modeltype, sequence, verbose = False):

	if (verbose):
		print(f"Making predictions using {modeltype}:")
		print(f"Sequence: {sequence}")

	encoding = utils.generate_data(sequence)['predict'][0]
	prediction = lmodel.predict(encoding)

	if ("PolyaClassifier" in modeltype):
		if (verbose):
			print("PolyaClassifier probability:", round(prediction[0][0],6))
		else:
			print(prediction[0][0])
		
	elif ("PolyaCleavage" in modeltype):
		if (verbose):
			print("PolyaCleavage distribution:", prediction[0].flatten().round(6).tolist())
		else:
			print(prediction[0].flatten())
		
	elif ("PolyaStrength" in modeltype):
		if (verbose):
			print("PolyaStrength score:", round(prediction[0][0],6))
		else:
			print(prediction[0][0])

	return



##############################
## DEFINE INPUTS
##############################

parser     = argparse.ArgumentParser(description = 'Make new predictions using yeast PolyaModels')
subparsers = parser.add_subparsers(help = 'sub-command help', dest = 'input_type')

parser_position = subparsers.add_parser('from_position', help = 'from position help')
parser_position.add_argument('-p', '--position',   type = str, required = True, help = 'String indicating a single genomic position for prediction. Sites must be formatted as chromosome:position:strand, where the position is the same as the start coordinate in a BED file.')
parser_position.add_argument('-g', '--genome',     type = str, required = True, help = 'Absolute path to a genome FASTA file containing the chromosome sequences.')
parser_position.add_argument('-c', '--chromsizes', type = str, required = True, help = 'Absolute path to a file containing the chromosome names and sizes in nucleotides. The file should be two columns.')
parser_position.add_argument('-m', '--model',      type = str, choices = ['S.cer PolyaClassifier','S.cer PolyaCleavage','S.cer PolyaStrength','S.pom PolyaClassifier'], help = 'Name of model to be used for predictions.')
parser_position.add_argument(      '--verbose',    action = 'store_true',       help = 'Flag indicating that verbose output should be printed.')

parser_sequence = subparsers.add_parser('from_sequence', help = 'from sequence help')
parser_sequence.add_argument('-s', '--sequence',   type = str, default = None, help = 'String containing 500 nt sequence to be directly used for prediction. This sequence should not contain any N nucleotides.')
parser_sequence.add_argument('-m', '--model',      type = str, choices = ['S.cer PolyaClassifier','S.cer PolyaCleavage','S.cer PolyaStrength','S.pom PolyaClassifier'], help = 'Name of model to be used for predictions.')
parser_sequence.add_argument(      '--verbose',    action = 'store_true',       help = 'Flag indicating that verbose output should be printed.')

args = parser.parse_args()



## PREPARE MODELS

if (args.model == 'S.cer PolyaClassifier'):
	lmodel = Ensemble([f"S_cerevisiae.PolyaClassifier.model_{i}.h5" for i in range(1,4)])

elif (args.model == 'S.cer PolyaCleavage'):
	lmodel = Predictor("S_cerevisiae.PolyaCleavage.h5")

elif (args.model == 'S.cer PolyaStrength'):
	lmodel = Predictor("S_cerevisiae.PolyaStrength.h5")

elif (args.model == 'S.pom PolyaClassifier'):
	lmodel = Ensemble([f"S_pombe.PolyaClassifier.model_{i}.h5" for i in range(1,4)])



## MAKE PREDICTIONS FOR INPUT SEQUENCE

if (args.input_type == 'from_sequence'):

	sequence = args.sequence.upper()

	if ('N' in sequence):
		raise ValueError("Input sequence contains missing nucleotides.")

	if (len(sequence) != 500):
		raise ValueError("Input sequence is not 500 nt.")

	report_predictions(lmodel, args.model, sequence, verbose = args.verbose)



## MAKE PREDICTIONS FOR INPUT POSITION

elif (args.input_type == 'from_position'):

	genome = utils.load_genome(args.genome)
	chrom_sizes = utils.get_chrom_size(args.chromsizes)

	chrom, position, strand = [int(x) if (i == 1) else str(x) for i,x in enumerate(args.position.split(":"))]
	sequence = utils.get_window(genome, chrom_sizes, chrom, position, strand)

	if ('N' in sequence):
		raise ValueError("Input sequence contains missing nucleotides.")

	if (len(sequence) != 500):
		raise ValueError("Input sequence is not 500 nt.")

	report_predictions(lmodel, args.model, sequence, verbose = args.verbose)




