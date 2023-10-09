#!/usr/bin/env python
# coding: utf-8

import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'



def load_genome(genome_path):
	'''Loads searchable and indexed genome from absolute path for genome FASTA file.
	'''

	import pyfaidx
	return pyfaidx.Fasta(genome_path, sequence_always_upper=True)



def get_chrom_size(chrom_size_path):
	'''Creates dictionary of chromosome sizes in nucleotides.
	'''

	chrom_dict = {}

	with open(chrom_size_path, mode = 'r') as infile:
		for line in infile:
			entries = line.strip("\n").split("\t")
			chrom_dict[entries[0]] = int(entries[1])

	return chrom_dict



def extract_sequence(genome, chrom, start, end, strand):
	'''Use genomic coordinates and indexed genome to extract the genomic sequence for the indicated interval.
	'''

	sequence = genome[chrom][start:end]

	if (strand == "+"):
		return sequence.seq.upper()
	elif (strand == '-'):
		return sequence.reverse.complement.seq.upper()



def get_window(genome, chrom_sizes, chrom, position, strand):
	'''Fetches 240 nt genomic sequence surrounding the indicated position.
	'''

	start = int(position - 120)
	end = int(position + 120)

	if (start <= 0) | (end >= chrom_sizes[chrom]):
		raise ValueError(f'Requested input with interval at ({chrom}:{start}-{end}:{strand}) exceeds the chromosome edges.')

	## if on the reverse strand, shift the coordinates one to put the position of interest at the 120th index in the sequence

	if (strand == '-'):
		start += 1
		end += 1

	return extract_sequence(genome, chrom, start, end, strand).upper()



def generate_data(sequences):
	'''Prepares data generator to flow data into the models for prediction.
	'''

	import numpy as np
	import pandas as pd
	import isolearn.keras as iso

	if (isinstance(sequences, str)):
		sequences = [sequences]

	df = pd.DataFrame.from_dict({'Sequence' : sequences})
	data_idx = np.arange(len(df), dtype = np.int)
	
	allSequenceGenerator = {
			gen_id : iso.DataGenerator(
				idx,
				{
					'df' : df
				},
				batch_size = len(idx),
				inputs = [
					{
						'id'          : 'seq',
						'source_type' : 'dataframe',
						'source'      : 'df',
						'extractor'   : lambda row,index: row['Sequence'].upper(),
						'encoder'     : iso.OneHotEncoder(seq_length = 240),
						'dim'         : (240,4),
						'sparsify'    : False,
					},
				],
				randomizers = [],
				shuffle = False,
				densify_batch_matrices = True,
				move_outputs_to_inputs = False
			)
			for gen_id, idx in [('predict', data_idx)]
		}
	return allSequenceGenerator


