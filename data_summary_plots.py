import json
import sys
import time
import matplotlib.pyplot as plt

t0 = time.perf_counter()
with open(sys.argv[1]) as fp:
	data = json.load(fp)
	t1 = time.perf_counter()
	print(f'time loading data: {t1-t0:.4f}')
fp.close()

print(len(data))

ids = list()
sizes = list()
shiftnum = list()
shiftp = list()
total = 0

for row in data:
	
	sizes.append(len(row['seq']))
	
	nshifts = 0
	for h, ha in zip(row['H'],row['HA']):
		if h is not None: nshifts += 1
		if ha is not None: nshifts += 1
	
	p = nshifts/(2*len(row['seq']))
	
	shiftnum.append(nshifts)
	shiftp.append(p)
	total += nshifts
	
	if row['name'] not in ids:
		ids.append(row['name'])

print(len(ids))
print(total)

plt.hist(sizes)
plt.show()

plt.hist(shiftnum)
plt.show()

plt.hist(shiftp)
plt.show()