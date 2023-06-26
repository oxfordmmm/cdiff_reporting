#!/usr/bin/env python3

import json, glob, sys, time, bisect, os
from optparse import OptionParser
import pandas as pd


def compareHashes(compareOut, hash_folder):
	#list of bad genes to exclude
	excludeList = ['CD630_17960', 'CD630_23930', 'CD630_25030', 'CD630_08260', 'CD630_09010', 'CD630_12750', 'CD630_22750', 'CD630_27680', 'CD630_11800', 'CD630_12080', 'CD630_15400', 'CD630_16950', 'CD630_20730', 'CD630_32240', 'CD630_34030']
	#fileList = glob.glob('%s/*.json'%input_folder)
	
	fileListhash = glob.glob('%s/*.json'%hash_folder)
	jsonList = []
	
	sys.stdout.write("Reading in JSON files\n")
	start = time.time()
	
	ids = []
	for f in fileListhash:
		with open(f, 'r') as fp:
			j = json.load(fp)
			id = os.path.basename(f).replace('_cgmlst.json', '')
			j['id'] = id
			ids.append(id)
			jsonList.append(j)
	end = time.time()
	sys.stdout.write("Seconds to read in files: %s\n"%(end - start))

	comparisons = pd.DataFrame(index = ids, columns = ids)
	diffs = pd.DataFrame(index = ids, columns = ids)
	distances = pd.DataFrame(index = ids, columns = ids)
	
	sys.stdout.write("Comparing profiles\n")
	start = time.time()

	for i in range(0, len(jsonList)):
		print(i)
		for j in range(i+1, len(jsonList)):
			a1 = [jsonList[i]['alleles'][k] for k in jsonList[i]['alleles'].keys() if k not in excludeList]
			a2 = [jsonList[j]['alleles'][k] for k in jsonList[j]['alleles'].keys() if k not in excludeList]
			id1 = jsonList[i]['id']
			id2 = jsonList[j]['id']
			compared = len([True for l1,l2 in zip(a1,a2) if l1 and l2])
			diff = len([True for l1,l2 in zip(a1,a2) if l1 and l2 and l1!=l2])
			
			distance = round(diff/compared, 3) if compared > 0 else 100000

			comparisons.loc[id1, id2] = compared
			diffs.loc[id1, id2] = diff
			distances.loc[id1, id2] = distance

			comparisons.loc[id2, id1] = compared
			diffs.loc[id2, id1] = diff
			distances.loc[id2, id1] = distance
	
	end = time.time()
	sys.stdout.write("Seconds to compare profiles: %s\n"%(end - start))

	comparisons.to_csv(f"{compareOut}_comps.csv")
	diffs.to_csv(f"{compareOut}_diffs.csv")
	distances.to_csv(f"{compareOut}_distances.csv")
	


if __name__ == "__main__":
	
	parser = OptionParser()
	parser.add_option( '-f', '--hash_folder', action = 'store', type='string', dest = 'hash_folder', default = '.' )
	parser.add_option( '-o', '--output_prefix', action = 'store', type='string', dest = 'output', default = 'output' )
	
	opts, args = parser.parse_args()
	
	compareHashes(opts.output, opts.hash_folder)
	