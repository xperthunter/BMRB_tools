#!/usr/bin/python3

import os
import sys
import re
import argparse
import json

import pynmrstar

parser = argparse.ArgumentParser(description='BMRB STAR3.0 to json')
parser.add_argument('--folder', required=True, type=str,
	metavar='<str>', help='File path to bmrb entries')
parser.add_argument('--json', required=True, type=str,
	metavar='<str>', help='Name of json file to dump to')
parser.add_argument('--err', required=True, type=str,
	metavar='<str>', help='file name for error output')
parser.add_argument('--csv', required=False, type=str,
	metavar='<str>', help='input csv file with specific ids to make json from')                

# grabbing all the args                                                          
arg = parser.parse_args()


aa_dict = {
	'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 
	'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L', 
	'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R',
	'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y'
}

aa_list = [
	'A', 'C', 'D', 'E', 'F', 
	'G', 'H', 'I', 'K', 'L', 
	'M', 'N', 'P', 'Q', 'R',
	'S', 'T', 'V', 'W', 'Y',
	'X'
]

class bmrbError(Exception):
	pass

def csv_to_list(file):
	with open(file, mode='r') as fp:
		ids = list()
		for line in fp.readlines():
			if re.search('[a-zA-Z]', line): continue
			
			ids.append(line.rstrip())
	
	return ids
		
def assembly_number(entry):
	
	assemblies = entry.get_tag('_Assembly.ID')
	
	return len(assemblies)
	
def get_shiftlist(entry):
	
	shift_loops = entry.get_loops_by_category("Atom_chem_shift")
	
	return shift_loops

def get_ents_with_shifts(shifts):
	eids = shifts.get_tag('Entity_ID')
	
	data_ents = []
	
	for ei in eids:
		if ei not in data_ents:
			data_ents.append(ei)
			continue
	
	if len(data_ents) == 0: return None
	
	return data_ents

def seq_by_entity(entry, data_ents):
	
	entity_sequences = dict()
	
	ent_info = entry.get_tags(['_Entity_assembly.Entity_ID', 
							   '_Entity_assembly.Entity_label'])
	
	found = dict()
	for id, label in zip(ent_info['_Entity_assembly.Entity_ID'],ent_info['_Entity_assembly.Entity_label']):
		if id not in data_ents:
			continue
		else:
			if id in found:
				continue
			else:
				found[id] = 1
		
		label = label.replace('$', '')
		entity_frame = entry.get_saveframe_by_name(label)
		
		seq = entity_seq(entity_frame)
		if seq is None: return None

		entity_sequences[id] = seq
	
	if len(found.keys()) == 0:
		print(data_ents)
		print(ent_info)
		print('never found entities with data')
		sys.exit()
		return None
	
	return entity_sequences

def entity_seq(entity_frame):
	polymer = entity_frame['Polymer_type']
	assert(len(polymer) == 1)
	
	if polymer[0] == '.':
		seq = entity_frame['Polymer_seq_one_letter_code']
		
		if len(seq) != 1:
			print(seq)
			return None
		elif seq[0] == '.':
			print(seq)
			return None
		else:
			seq = seq[0].replace('\n','')
			seq = seq.rstrip()
			seq = seq.upper()
			alphabet = dict()
			for aa in seq:
				if aa in alphabet: continue
				else: alphabet[aa] = 1
				
			if len(alphabet.keys() > 5):
				return seq[0].replace('\n', '')
			else:
				print(seq)
				return None
	
	elif polymer[0] == 'polypeptide(L)':
		seq = entity_frame['Polymer_seq_one_letter_code']
		
		if len(seq) != 1:
			print(seq)
			return None
		
		seq = seq[0].replace('\n','')
		seq = seq.rstrip()
		seq = seq.upper()
		return seq
	else:
		print(polymer[0])
		return None

def get_shifts(shifts, proteins, id):
	eids = shifts.get_tag('Entity_ID')
	ind = shifts.get_tag('Seq_ID')
	aas = shifts.get_tag('Comp_ID')
	atm = shifts.get_tag('Atom_ID')
	val  = shifts.get_tag('Val')
	
	assert(len(eids) == len(ind) == len(aas) == len(atm) == len(val))
	
	shift_data = dict()
	for id in proteins.keys():
		shift_data[id] = dict()
		shift_data[id]['HA'] = []
		shift_data[id]['H'] = []
		HA = []
		H = []
		for i in range(len(proteins[id])):
			HA.append(None)
			H.append(None)
		
		shift_data[id]['HA'] = HA
		shift_data[id]['H'] = H
	
	for i in range(len(eids)):
		id = eids[i]
		pos = int(ind[i])-1
		if aas[i] not in aa_dict:
			res = 'X'
		else:
			res = aa_dict[aas[i]]
		atom = atm[i]
		cs = float(val[i])
		
		if pos >= len(proteins[id]):
#			print(f'{pos} {proteins[id]}')
			return None
		
		if proteins[id][pos] != res:
#			print(f'{pos} {res}')
			return None
		
		if atom == 'HA':
			shift_data[id]['HA'][pos] = cs
		elif atom == 'H':
			shift_data[id]['H'][pos] = cs
	
	return shift_data

weird = {
	'STAR' : [],
	'assembly>1' : [],
	'!shiftloops' : [],
	'shiftlist!1' : [],
	'seqsError' : [],
	'shifterrors' : [],
	'ents_with_cs_error' : [],
	'not_requested' : []
}

if arg.csv:
	id_set = csv_to_list(arg.csv)
#	print(id_set)
#	sys.exit()

data = list()

for file in os.listdir(arg.folder):
	if not re.search(r'str$', file): continue
	try:
		ent = pynmrstar.Entry.from_file(f'{arg.folder}/{file}')
	except:
		weird['STAR'].append(file)
		continue
	print(file)
	
	id = ent.get_tag('_Entry.ID')
	bid = id[0]
	
	if arg.csv:
		if bid not in id_set:
			weird['not_requested'].append(bid)
			continue
	
	if assembly_number(ent) != 1:
		weird['assembly>1'].append(file)
		continue
	
	shift_sets = get_shiftlist(ent)

	if len(shift_sets) == 0:
		weird['!shiftloops'].append(file)
		continue
	
	if len(shift_sets) != 1:
		weird['shiftlist!1'].append(file)
		continue
	
	entities_with_cs = get_ents_with_shifts(shift_sets[0])
	if len(entities_with_cs) is None:
		weird['ents_with_cs_error'].append(file)
		continue
	
	seqs_per_ent = seq_by_entity(ent,entities_with_cs)
	if seqs_per_ent is None:
		weird['seqsError'].append(file)
		continue
	
	cs_data = get_shifts(shift_sets[0], seqs_per_ent, bid)
	
	if cs_data is None:
		weird['shifterrors'].append(file)
		continue
	
	for id in cs_data:
		data.append({'name' : bid, 'seq' : seqs_per_ent[id], 'H' : cs_data[id]['H'], 'HA' : cs_data[id]['HA']})
	
	print('good')
	
with open(arg.json, 'w') as fp:
	json.dump(data, fp)
fp.close()

with open(arg.err, 'w') as fe:
	json.dump(weird, fe)
fe.close()	


"""
def one_component(entry):
	(num_comp,) = entry.get_tag('_Assembly.Number_of_components')
	
	if num_comp == '1':
		return True
	else:
		entities = entry.get_tag('_Entity_assembly.Entity_label')
		if len(entities) == 1:
			return True
		else:
			return False

def get_protein(entry):
	polymers = entry.get_tag('_Entity.Polymer_type')
	
	if polymers[0] != 'polypeptide(L)':
		return None
	
	seq = entry.get_tag('_Entity.Polymer_seq_one_letter_code')
	
	return seq[0].replace('\n','')

def get_seq(entry):
	
	entities = entry.get_tags(['_Entity_assembly.Entity_label',
							  '_Entity_assembly.Experimental_data_reported'])
	
	counter = 1
	for sframe in entities['_Entity_assembly.Entity_label']:
		sframe = sframe.replace('$', '')
		frame = entry.get_saveframe_by_name(sframe)
		type = frame['Type']
		polymer_type = frame['Polymer_type']
		seq = frame['Polymer_seq_one_letter_code']
		if seq[0] == ".":
			seq = frame['Nonpolymer_comp_ID'][0]
		else:
			seq = seq[0].replace('\n','')
		
		yield counter, type[0], polymer_type[0], seq, entities['_Entity_assembly.Experimental_data_reported'][counter-1]
		counter += 1

	return
"""
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	