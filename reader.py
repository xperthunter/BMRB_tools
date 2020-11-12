import json
import sys
import time
import lzma
import pandas as pd

t0 = time.perf_counter()
df = pd.read_json(sys.argv[1]).reset_index()
t1 = time.perf_counter()
print(f'time loading data: {t1-t0:.4f}')

print(df.head(5))
print(df.columns.values)
"""
t0 = time.perf_counter()
with lzma.open(sys.argv[1]) as fp:
	data = json.load(fp)
	t1 = time.perf_counter()
	print(f'time loading data: {t1-t0:.4f}')
fp.close()

print(len(data))
"""
t0 = time.perf_counter()
with open(sys.argv[2]) as fe:
	err = json.load(fe)
	t1 = time.perf_counter()
	print(f'time loading errors: {t1-t0:.4f}')
fe.close()

for k in err.keys():
	print(k, len(err[k]))

print(df.isnull().sum().sum())

#for index, row in df.iterrows():
#	print(row)
#	print (row.isnull().any())
#	for c in df.columns.values:
#		print(row[c])
#		if row[c] == 'NaN': 
#			print(c)
#			sys.exit()
#	sys.exit()
"""
	rsize = round(len(row['seq']),-1)
	if rsize in sizes:
		sizes[rsize] += 1
	else:
		sizes[rsize] = 1
	
	nshifts = 0
	for h, ha in zip(row['H'],row['HA']):
		if h is not None: nshifts += 1
		if ha is not None: nshifts += 1
	
	p = nshifts/(2*len(row['seq']))
	p = round(p, 1)
	
	if p in shiftp:
		shiftp[p] += 1
	else:
		shiftp[p] = 1
	
	nshifts = round(nshifts, -1)
	if nshifts in shiftnum:
		shiftnum[nshifts] += 1
	else:
		shiftnum[nshifts] = 1 
	
	if row['name'] not in ids:
		ids.append(row['name'])

print(len(ids))

for s in sorted(sizes.keys()):
	print(s, sizes[s])
print()
for n in sorted(shiftnum.keys()):
	print(n, shiftnum[n])
print()
for f in sorted(shiftp.keys()):
	print(f, shiftp[f])
"""