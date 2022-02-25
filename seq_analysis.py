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
"""
parser.add_argument('--json', required=True, type=str,
	metavar='<str>', help='Name of json file to dump to')
parser.add_argument('--err', required=True, type=str,
	metavar='<str>', help='file name for error output')                          
"""

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

protein_only = [
	'D', 'E', 'F', 'H', 'I',
	'K', 'L', 'M', 'N', 'P',
	'Q', 'R', 'S', 'V', 'W',
	'Y'
]

class bmrbError(Exception):
	pass

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
	
	if polymer[0] == 'polypeptide(L)':
		seq = entity_frame['Polymer_seq_one_letter_code']
		
		if len(seq) != 1: return None
		
		if seq[0] == '.':
			seq = entity_frame['Polymer_seq_one_letter_code_can']
			if len(seq) != 1: return None
			
			if seq[0] == '.': return None
		
		seq = seq[0].replace('\n','')
		seq = seq.rstrip()
		seq = seq.upper()
				
		return seq
	
	elif polymer[0] == '.':
		seq = entity_frame['Polymer_seq_one_letter_code']
		if len(seq) != 1: return None
			
		if seq[0] == '.':
			seq = entity_frame['Polymer_seq_one_letter_code_can']
			if len(seq) != 1: return None
			
			if seq[0] == '.': return None
			
		seq = seq[0].replace('\n','')
		seq = seq.rstrip()
		seq = seq.upper()
		for aa in seq:
			if aa in protein_only:
				return seq
		
		return None
		
	else:
		print(polymer[0])
		return None

def seq_check(seq):
	print(f'checking {seq}')
	for letter in seq:
		if letter == ' ': 
			print(letter, seq)
			sys.exit()
		if letter in aa_list: continue
		else:
			print(letter, seq)
			sys.exit()
	return True

weird = {
	'STAR' : [],
	'assembly>1' : [],
	'!shiftloops' : [],
	'shiftlist!1' : [],
	'seqsError' : [],
	'shifterrors' : [],
	'ents_with_cs_error' : []
}

data = list()

for file in os.listdir(arg.folder):
	if not re.search(r'str$', file): continue
	path = f'{arg.folder}{file}'
	if os.path.islink(path):
		print(os.readlink(path))
		path = os.path.join(os.path.dirname(arg.folder), os.readlink(path))
	try:
		ent = pynmrstar.Entry.from_file(path)
	except:
		weird['STAR'].append(file)
		continue
	print(file)
	
	id = ent.get_tag('_Entry.ID')
	bid = id[0]
	
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
	
	for s in seqs_per_ent.values():
		if seq_check(s): 
			print(f'>{s}*/*/*/*')
			continue
		else: sys.exit()
	
	"""
	cs_data = get_shifts(shift_sets[0], seqs_per_ent, bid)
	
	if cs_data is None:
		weird['shifterrors'].append(file)
		continue
	
	for id in cs_data:
		data.append({'name' : bid, 'seq' : seqs_per_ent[id], 'H' : cs_data[id]['H'], 'HA' : cs_data[id]['HA']})
	
	print('good')
	"""

print(json.dumps(weird, indent=4))