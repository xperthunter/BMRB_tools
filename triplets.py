
import json
import sys
import statistics

data = None
with open('refdb.json') as fp:
	data = json.load(fp)

count = {}
for item in data:
	seq = item['seq']
	for i in range(1, len(seq) -1):
		if seq[i] != 'L': continue # focus on most common first
		a3 = seq[i-1:i+2]
		
		for atom in item['cs'][i]:
			if item['cs'][i][atom] != None:
				if atom not in count: count[atom] = {}
				if a3 not in count[atom]: count[atom][a3] = []
				count[atom][a3].append(item['cs'][i][atom])

for atom in count:
	filename = 'cs_' + atom + '.txt'
	with open(filename, "w") as fp:
		for a3 in count[atom]:
			mean = statistics.mean(count[atom][a3])
			fp.write(f'{atom}\t{a3}\t{len(count[atom][a3])}\t{mean:.3f}\n')

