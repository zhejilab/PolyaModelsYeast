#!/usr/bin/env python
# coding: utf-8


## IMPORTS

import yeast_utilities as utils

import os, argparse
import numpy as np

from contextlib import redirect_stderr

with redirect_stderr(open(os.devnull, "w")):
	from keras.models import load_model



## DEFINE INPUTS

parser     = argparse.ArgumentParser(description = 'Make new predictions using yeast PolyaModels')
subparsers = parser.add_subparsers(help = 'sub-command help', dest = 'input_type')

parser_position = subparsers.add_parser('from_position', help = 'from position help')
parser_position.add_argument('-p', '--position',   type = str, required = True, help = 'String indicating a single genomic position for prediction. Sites must be formatted as chromosome:position:strand, where the position is the same as the start coordinate in a BED file.')
parser_position.add_argument('-g', '--genome',     type = str, required = True, help = 'Absolute path to a genome FASTA file containing the chromosome sequences.')
parser_position.add_argument('-c', '--chromsizes', type = str, required = True, help = 'Absolute path to a file containing the chromosome names and sizes in nucleotides. The file should be two columns.')
parser_position.add_argument('-m', '--model',      type = str, choices = ['S.cer PolyaClassifier','S.cer PolyaCleavage','S.cer PolyaStrength','S.pom PolyaClassifier'], help = 'Name of model to be used for predictions.')

parser_sequence = subparsers.add_parser('from_sequence', help = 'from sequence help')
parser_sequence.add_argument('-s', '--sequence',   type = str, default = None, help = 'String containing 240 nt sequence to be directly used for prediction. This sequence should not contain any N nucleotides.')
parser_sequence.add_argument('-m', '--model',      type = str, choices = ['S.cer PolyaClassifier','S.cer PolyaCleavage','S.cer PolyaStrength','S.pom PolyaClassifier'], help = 'Name of model to be used for predictions.')

args = parser.parse_args()



## PREPARE MODELS

if (args.model == 'S.cer PolyaClassifier'):
	lmodel = load_model("S_cerevisiae.PolyaClassifier.h5")
elif (args.model == 'S.cer PolyaCleavage'):
	lmodel = load_model("S_cerevisiae.PolyaCleavage.h5")
elif (args.model == 'S.cer PolyaStrength'):
	lmodel = load_model("S_cerevisiae.PolyaStrength.h5")
elif (args.model == 'S.pom PolyaClassifier'):
	lmodel = load_model("S_pombe.PolyaClassifier.h5")



## MAKE PREDICTIONS FOR INPUT SEQUENCE

if (args.input_type == 'from_sequence'):

	sequence = args.sequence.upper()

	if ('N' in sequence):
		raise ValueError("Input sequence contains missing nucleotides.")

	if (len(sequence) != 240):
		raise ValueError("Input sequence is not 240 nt.")

	print("Sequence:", sequence)

	encoding = utils.generate_data(sequence)['predict'][0]
	prediction = lmodel.predict(encoding)

	if ("PolyaClassifier" in args.model):
		print("PolyaClassifier probability:", round(prediction[0][0],6))
		
	elif ("PolyaCleavage" in args.model):

		rawcleavage = prediction[0].flatten()

		subtracted = rawcleavage - 0.02
		subtracted[subtracted <= 0] = 0
		normcleavage = subtracted / np.sum(subtracted) if (np.sum(subtracted) > 0) else np.asarray([0]*50)

		print("PolyaCleavage distribution:", normcleavage.round(6).tolist())
		
	elif ("PolyaStrength" in args.model):
		print("PolyaStrength score:", round(prediction[0][0],6))



## MAKE PREDICTIONS FOR INPUT POSITION

elif (args.input_type == 'from_position'):

	genome = utils.load_genome(args.genome)
	chrom_sizes = utils.get_chrom_size(args.chromsizes)

	chrom, position, strand = [int(x) if (i == 1) else str(x) for i,x in enumerate(args.position.split(":"))]
	sequence = utils.get_window(genome, chrom_sizes, chrom, position, strand)

	if ('N' in sequence):
		raise ValueError("Input sequence contains missing nucleotides.")

	if (len(sequence) != 240):
		raise ValueError("Input sequence is not 240 nt.")

	print("Sequence:", sequence)

	encoding = utils.generate_data(sequence)['predict'][0]
	prediction = lmodel.predict(encoding)

	if ("PolyaClassifier" in args.model):
		print("PolyaClassifier probability:", round(prediction[0][0],6))
		
	elif ("PolyaCleavage" in args.model):

		rawcleavage = prediction[0].flatten()

		subtracted = rawcleavage - 0.02
		subtracted[subtracted <= 0] = 0
		normcleavage = subtracted / np.sum(subtracted) if (np.sum(subtracted) > 0) else np.asarray([0]*50)

		print("PolyaCleavage distribution:", normcleavage.round(6).tolist())
		
	elif ("PolyaStrength" in args.model):
		print("PolyaStrength score:", round(prediction[0][0],6))


