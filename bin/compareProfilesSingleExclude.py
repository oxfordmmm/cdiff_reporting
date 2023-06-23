#!/usr/bin/env python3

import json, glob, sys, time, bisect
from optparse import OptionParser


def compareHash(input_hash, compareOut, hash_folder, distance_cutoff, min_distance_cutoff, min_number_of_samples):
	#list of bad genes to exclude
	excludeList = ['CD630_17960', 'CD630_23930', 'CD630_25030', 'CD630_08260', 'CD630_09010', 'CD630_12750', 'CD630_22750', 'CD630_27680', 'CD630_11800', 'CD630_12080', 'CD630_15400', 'CD630_16950', 'CD630_20730', 'CD630_32240', 'CD630_34030']
	#fileList = glob.glob('%s/*.json'%input_folder)
	
	fileListhash = glob.glob('%s/*.json'%hash_folder)
	jsonList = []
	
	sys.stdout.write("Reading in JSON files\n")
	start = time.time()

	in_hash = None
	with open(input_hash, 'r') as fp:
		in_hash = json.load(fp)
	
	for f in fileListhash:
		with open(f, 'r') as fp:
			j = json.load(fp)
			jsonList.append(j)
	end = time.time()
	sys.stdout.write("Seconds to read in files: %s\n"%(end - start))
	
	closest_samples = []
	sys.stdout.write("Comparing profiles\n")
	start = time.time()

	sample_comparisons = []
	for i in range(0, len(jsonList)):
		a1 = [jsonList[i]['alleles'][k] for k in jsonList[i]['alleles'].keys() if k not in excludeList]
		a2 = [in_hash['alleles'][k] for k in in_hash['alleles'].keys() if k not in excludeList]
		compared = len([True for l1,l2 in zip(a1,a2) if l1 and l2])
		diff = len([True for l1,l2 in zip(a1,a2) if l1 and l2 and l1!=l2])
		
		if compared > 0:
			sample = (jsonList[i]['name'], compared, diff, "%0.3f"%(diff/compared))
			sample_comparisons.append((sample, diff))
	
	end = time.time()
	sys.stdout.write("Seconds to compare profiles: %s\n"%(end - start))
	start = time.time()

	# sort by differences ascending
	sample_comparisons = sorted(sample_comparisons, key=lambda x: x[1])
	# Filter out ones which are too close
	sample_comparisons = [(sample, diff) for sample, diff in sample_comparisons if diff > min_distance_cutoff]
	# count number which are close enough
	num_close = len([1 for sample, diff in sample_comparisons if diff < distance_cutoff])
	sample_comparisons = sample_comparisons[:max(min_number_of_samples, num_close)]
	closest_samples = [sample for sample, diff in sample_comparisons]

	with open(compareOut, 'w') as w:
		w.write('Sample\tLoci_compared\tDifferences\tDistance\n')
		for sample in closest_samples:
			w.write('%s\t%s\t%s\t%s\n'%(sample[0], sample[1], sample[2], sample[3]))
	
	end = time.time()
	sys.stdout.write("Seconds to sort: %s\n"%(end - start))

if __name__ == "__main__":
	
	parser = OptionParser()
	parser.add_option( '-f', '--hash_folder', action = 'store', type='string', dest = 'hash_folder', default = '.' )
	parser.add_option( '-i', '--input_hash', action = 'store', type='string', dest = 'input_hash', default = '.' )
	parser.add_option( '-d', '--distance_cutoff', action = 'store', type='int', dest = 'distance_cutoff', default = 20)
	parser.add_option( '-m', '--min_distance_cutoff', action = 'store', type='int', dest = 'min_distance_cutoff', default = 1)
	parser.add_option( '-n', '--min_number_of_samples', action = 'store', type='int', dest = 'min_number_of_samples', default = -1,
							help = "Number of samples to return regardless of distance")
	parser.add_option( '-o', '--output', action = 'store', type='string', dest = 'output', default = 'output' )
	
	opts, args = parser.parse_args()
	
	compareHash(opts.input_hash, opts.output, opts.hash_folder, opts.distance_cutoff, opts.min_distance_cutoff, opts.min_number_of_samples)
	
#E.g. compareProfiles.py -i /well/bag/deyre/analysis/spades-flow/replicates_output -o  /well/bag/deyre/analysis/spades-flow/replicates_compare.txt

# compareProfiles.py -i /home/davideyre/hash-cgmlst/comparison_study_data/replicates_output -o  /home/davideyre/hash-cgmlst/comparison_study_data/replicates_compare.txt

