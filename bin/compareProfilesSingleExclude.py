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
	for i in range(0, len(jsonList)):
		a1 = [jsonList[i]['alleles'][k] for k in jsonList[i]['alleles'].keys() if k not in excludeList]
		a2 = [in_hash['alleles'][k] for k in in_hash['alleles'].keys() if k not in excludeList]
		compared = 0
		diff = 0
		inserted = False
		compared = len([True for l1,l2 in zip(a1,a2) if l1 and l2])
		diff = len([True for l1,l2 in zip(a1,a2) if l1 and l2 and l1!=l2])
		
		sample = (jsonList[i]['name'], compared, diff, "%0.3f"%(diff/compared))
		if diff < distance_cutoff:
			if diff > min_distance_cutoff:
				if len(closest_samples) > 0:
					closest_samples.insert(bisect.bisect_left([i[2] for i in closest_samples], sample[2]), sample)
				else:
					closest_samples = [sample]
				inserted = True
		elif min_number_of_samples > 0:
			if len(closest_samples) > 0 and len(closest_samples) < min_number_of_samples:
				index = bisect.bisect_left([i[2] for i in closest_samples], sample[2])
				if index < min_number_of_samples:
					closest_samples.insert(index, sample)
				else:
					closest_samples = [sample]
				inserted = True
		
		if inserted and len(closest_samples) > min_number_of_samples:
			while len(closest_samples) > min_number_of_samples:
				if closest_samples[-1][2] < distance_cutoff:
					break
				closest_samples.pop()


	with open(compareOut, 'w') as w:
		w.write('Sample\tLoci_compared\tDifferences\tDistance\n')
		for sample in closest_samples:
			w.write('%s\t%s\t%s\t%s\n'%(sample[0], sample[1], sample[2], sample[3]))
	
	end = time.time()
	sys.stdout.write("Seconds to compare profiles: %s\n"%(end - start))


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

