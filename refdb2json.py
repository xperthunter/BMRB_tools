
import json
import sys
import re


def myfloat(s):
	if s == '****': return None
	else: return float(s)

def parse(line, type):
	m = re.search(': (.+)', line)
	if m:
		if   type == 's': return m[1]
		elif type == 'f': return float(m[1])
		elif type == 'i': return int(m[1])
	else: return None

def proc(desc, lines, type):
	obj = {
		'bid':  parse(lines[2], 'i'), # BMRB ID
		'pid':  parse(lines[1], 's'), # PDB ID
		'lit':  parse(lines[0], 's'), # literature reference
		'fold': parse(lines[3], 's'), # single letter ID
		'ph':   parse(lines[4], 'f'), # not always listed
		'temp': parse(lines[5], 'f'), # not always listed
		'slen': parse(lines[6], 'i'), # number of residues
		'seq':  None, # to be completed later
		'type': type, # N, C, H
		'data': []
	}
	
	for line in lines[7:]:
		f = line.split()
		if len(f) == 0: continue
		elif type == 'h':
			if len(f) != 13:
				obj['error'] = True
			else:
				obj['data'].append({
					'idx': f[0],
					'aa': f[1],
					'st': f[2],            # Trp example
					'v1': myfloat(f[3]),   # H
					'v2': myfloat(f[4]),   # HA
					'v3': myfloat(f[5]),   # HB2
					'v4': myfloat(f[6]),   # HB3
					'v5': myfloat(f[7]),   # HD1
					'v6': myfloat(f[8]),   # HE1
					'v7': myfloat(f[9]),   # HE3
					'v8': myfloat(f[10]),  # HZ2
					'v9': myfloat(f[11]),  # HZ3
					'v10': myfloat(f[12]), # HH2
				})
		elif type == 'c':
			if len(f) != 7:
				obj['error'] = True
			else:
				obj['data'].append({
					'idx': f[0],
					'aa': f[1],
					'st': f[2],
					'v1': myfloat(f[3]), # ? may be something with methionine?
					'v2': myfloat(f[4]), # CA
					'v3': myfloat(f[5]), # CB
					'v4': myfloat(f[6]), # C
				})
		elif type == 'n':
			if len(f) != 5:
				obj['error'] = True
			else:
				obj['data'].append({
					'idx': f[0],
					'aa': f[1],
					'st': f[2],
					'v1': myfloat(f[3]), # N
					'v2': myfloat(f[4]), # H
				})
		else:
			sys.stderr.write('unknown file type\n')
			sys.exit(0)
		
	return obj

def read_shifty(filename):
	desc = None
	lines = []
	type = None
	if   filename.endswith('C.db'): type = 'c'
	elif filename.endswith('H.db'): type = 'h'
	elif filename.endswith('N.db'): type = 'n'
	else:
		sys.stderr.write('unknown file inference\n')
		sys.exit(1)
	
	fp = None
	if filename == '-':
		fp = sys.stdin
	elif filename.endswith('.gz'):
		fp = gzip.open(filename, 'rt')
	else:
		fp = open(filename)

	for line in fp.readlines():
		if line.startswith('>'):
			if len(lines) > 0:
				yield(proc(desc, lines, type))
				desc = line[1:]
				lines = []
			else:
				desc = line[1:]
		else:
			lines.append(line)
	yield(proc(desc, lines, type))
	fp.close()

def get(stuff, field):
	found = {}
	for type in stuff:
		if field in stuff[type]:
			found[field] = stuff[type][field]
	if len(found) == 1:
		return found[field]
	else:
		sys.stderr.write(f'unexpected data mismatch\n')
		sys.exit(0)

def getseq(stuff):
	field = list(stuff.keys())
	seq = []
	sseq = []
	for d in stuff[field[0]]['data']:
		seq.append(d['aa'])
	return ''.join(seq)

def getstr(stuff):
	field = list(stuff.keys())
	seq = []
	for d in stuff[field[0]]['data']:
		seq.append(d['st'])
	return ''.join(seq)

def translate_nshift(data):
	cs = {}
	aa = data['aa']
	cs['N'] = data['v1']
	if   aa == 'R': cs['NE'] = data['v2']
	elif aa == 'K': cs['NZ'] = data['v2']
	return cs

def translate_hshift(data):
	cs = {}
	aa = data['aa']
	if aa == 'A':
		cs['H']  = data['v1']
		cs['HA'] = data['v2']
		cs['HB'] = data['v3']
	elif aa == 'C' or aa == 'D' or aa == 'E' or aa == 'S':
		cs['H']   = data['v1']
		cs['HA']  = data['v2']
		cs['HB2'] = data['v3']
		cs['HB3'] = data['v4']
	elif aa == 'F' or aa == 'Y':
		cs['H']   = data['v1']
		cs['HA']  = data['v2']
		cs['HB2'] = data['v3']
		cs['HB3'] = data['v4']
		cs['HD1'] = data['v7']
		cs['HD2'] = data['v8']
		cs['HE1'] = data['v9']
		cs['HE2'] = data['v10']
		#cs['HZ'] not recorded in SHIFTY file but is in RefDB
	elif aa == 'G':
		cs['H']   = data['v1']
		cs['HA2'] = data['v2']
		cs['HA3'] = data['v3']
	elif aa == 'H':
		cs['H']   = data['v1']
		cs['HA']  = data['v2']
		cs['HB2'] = data['v3']
		cs['HB3'] = data['v4']
		cs['HD2'] = data['v8']
		cs['HE1'] = data['v9']
	elif aa == 'I':
		cs['H']    = data['v1']
		cs['HA']   = data['v2']
		cs['HB']   = data['v3']
		cs['HG12'] = data['v4']
		cs['HG13'] = data['v5']
		cs['HG2']  = data['v6']
		cs['HD1']  = data['v7']
	elif aa == 'K':
		cs['H']   = data['v1']
		cs['HA']  = data['v2']
		cs['HB2'] = data['v3']
		cs['HB3'] = data['v4']
		cs['HG2'] = data['v5']
		cs['HG3'] = data['v6']
		cs['HD2'] = data['v7']
		cs['HD3'] = data['v8']
		cs['HE2'] = data['v9']
		cs['HE3'] = data['v10']
	elif aa == 'L':
		cs['H']   = data['v1']
		cs['HA']  = data['v2']
		cs['HB2'] = data['v3']
		cs['HB3'] = data['v4']
		cs['HG']  = data['v5']
		cs['HD1'] = data['v6']
		cs['HD2'] = data['v7']
	elif aa == 'M':
		cs['H']   = data['v1']
		cs['HA']  = data['v2']
		cs['HB2'] = data['v3']
		cs['HB3'] = data['v4']
		cs['HG2'] = data['v5']
		cs['HG3'] = data['v6']
		cs['HE']  = data['v7']
	elif aa == 'N':
		cs['H']    = data['v1']
		cs['HA']   = data['v2']
		cs['HB2']  = data['v3']
		cs['HB3']  = data['v4']
		cs['HD21'] = data['v7']
		cs['HD22'] = data['v8']
	elif aa == 'P':
		cs['HA']  = data['v2']
		cs['HB2'] = data['v3']
		cs['HB3'] = data['v4']
		cs['HG2'] = data['v5']
		cs['HG3'] = data['v6']
		cs['HD2'] = data['v7']
		cs['HD3'] = data['v8']
	elif aa == 'Q': 
		cs['H']    = data['v1']
		cs['HA']   = data['v2']
		cs['HB2']  = data['v3']
		cs['HB3']  = data['v4']
		cs['HG2']  = data['v5']
		cs['HG3']  = data['v6']
		cs['HE21'] = data['v9']
		cs['HE22'] = data['v10']
	elif aa == 'R':
		cs['H']   = data['v1']
		cs['HA']  = data['v2']
		cs['HB2'] = data['v3']
		cs['HB3'] = data['v4']
		cs['HG2'] = data['v5']
		cs['HG3'] = data['v6']
		cs['HD2'] = data['v7']
		cs['HD3'] = data['v8']
		cs['HE']  = data['v9']
	elif aa == 'T':
		cs['H']   = data['v1']
		cs['HA']  = data['v2']
		cs['HB']  = data['v3']
		cs['HG2'] = data['v5']
	elif aa == 'V':
		cs['H']   = data['v1']
		cs['HA']  = data['v2']
		cs['HB']  = data['v3']
		cs['HG1'] = data['v4']
		cs['HG2'] = data['v5']
	elif aa == 'W':
		cs['H']   = data['v1']
		cs['HA']  = data['v2']
		cs['HB2'] = data['v3']
		cs['HB3'] = data['v4']
		cs['HD1'] = data['v5']
		cs['HE1'] = data['v6']
		cs['HE3'] = data['v7']
		cs['HZ2'] = data['v8']
		cs['HZ3'] = data['v9']
		cs['HH2'] = data['v10']
	else:
		sys.stderr.write(f'unknown amino acid {aa}\n')
		cs['X'] = 0.0

	return cs

def translate_cshift(data):
	cs = {}
	aa = data['aa']
	
	if aa == 'G':
		cs['CA'] = data['v2']
	else: # ACDEFHIKLMNPQRSTVWY
		cs['C']  = data['v4']
		cs['CA'] = data['v2']
		cs['CB'] = data['v3']

	return cs

def getshifts(stuff):
	
	hdata = []
	if 'h' in stuff:
		for d in stuff['h']['data']:
			hdata.append(translate_hshift(d))

	cdata = []
	if 'c' in stuff:
		for d in stuff['c']['data']:
			cdata.append(translate_cshift(d))
	
	ndata = []
	if 'n' in stuff:
		for d in stuff['n']['data']:
			ndata.append(translate_nshift(d))	
	
	size = max(len(hdata), len(cdata), len(ndata))
	data = [ {} for i in range(size) ]
	for i in range(len(hdata)):
		for cs in hdata[i]:
			data[i][cs] = hdata[i][cs]
	for i in range(len(cdata)):
		for cs in cdata[i]:
			data[i][cs] = cdata[i][cs]
	for i in range(len(ndata)):
		for cs in ndata[i]:
			data[i][cs] = ndata[i][cs]
			
	return data
	

files = ('RefDB-C.db', 'RefDB-H.db', 'RefDB-N.db')
refdb = {}
for file in files:
	for entry in read_shifty(file):
		bid = entry['bid']
		type = entry['type']
		if bid not in refdb: refdb[bid] = {}
		if type not in refdb[bid]: refdb[bid][type] = entry

jdb = []
perrors = []
for bid in refdb:
	
	# check for errors
	error = False
	for type in refdb[bid]:
		if 'error' in refdb[bid][type]: error = True
	if error:
		perrors.append(bid)
		continue # skipping for now

	# build object
	obj = {
		'bmrb': bid,
		'pdb':  get(refdb[bid], 'pid'),
		'ref':  get(refdb[bid], 'lit'),
		'temp': get(refdb[bid], 'temp'),
		'ph':   get(refdb[bid], 'ph'),
		'seq':  getseq(refdb[bid]),
		'str':  getstr(refdb[bid]),
		'cs':   getshifts(refdb[bid]),
	}
	
	jdb.append(obj)

print(json.dumps(jdb, indent=4))
	

# error report
sys.stderr.write(f'parsing errors: {len(perrors)}\n')
sys.stderr.write(f'bmrb identifiers: {perrors}\n')





