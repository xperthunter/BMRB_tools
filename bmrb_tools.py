#!/usr/bin/env python3

import os
import sys
import re
import json

import pynmrstar

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


def csv_to_list(file):
	with open(file, mode='r') as fp:
		ids = list()
		for line in fp.readlines():
			if re.search('[a-zA-Z]', line): continue
			
			ids.append(line.rstrip())
	
	return ids


def seq_check(seq):
	for letter in seq:
		if letter == ' ': 
			print(letter, seq)
			sys.exit()
		if letter in aa_list: continue
		else:
			print(f'{letter} not in aa_list\n{seq}')
			sys.exit()
	return True


def assembly_number(entry):
	
	assemblies = entry.get_tag('_Assembly.ID')
	
	return len(assemblies)


def get_shiftlist(entry):
	
	shift_loops = entry.get_loops_by_category("Atom_chem_shift")
	
	return shift_loops


def conditions(entry):
	shiftframe = entry.get_saveframes_by_category("assigned_chemical_shifts")[0]
	condition_label = shiftframe['Sample_condition_list_label']
	if len(condition_label) != 1: return None
	else:                         condition_label = condition_label[0]
	condition_label = re.sub('^\$','',condition_label)
	
	cond_frame = entry.get_saveframe_by_name(condition_label)
	cond_loop = cond_frame.get_loop_by_category("sample_condition_variable")
	
	type  = cond_loop.get_tag('Type')
	val   = cond_loop.get_tag('Val')
	units = cond_loop.get_tag('Val_units')
	
	if len(type) == 0 or len(val) == 0 or len(units) == 0: return None
	
	conditions_dict = dict()
	for t,v,u in zip(type, val, units):
		conditions_dict[t] = {'val':v, 'units':u}
	
	return conditions_dict


def conformers_check(loop):
	model_id  = loop.get_tag('Model_ID')
	entity_id = loop.get_tag('Label_entity_ID')
	resid     = loop.get_tag('Label_comp_ID')
	
	max_model = int(model_id[-1])
	entsm = dict()
	seqsm = dict()
	entc  = dict()
	seqc  = dict()
	
	for m, e, r in zip(model_id, entity_id, resid):
		if m not in entsm:
			entsm[m] = dict()
			seqsm[m] = dict()
		
		if e not in entsm[m]:
			entsm[m][e] = True
			seqsm[m][e] = []
			if e not in entc:
				entc[e] = 0
			entc[e] += 1
					
		seqsm[m][e].append(r)
	
	for k, v in entc.items():
		if v != max_model: return False
	
	for k1 in seqsm.keys():
		for k2 in seqsm[k1].keys():
			seq = ''.join(seqsm[k1][k2])
			if seq not in seqc:
				seqc[seq] = 0
			seqc[seq] += 1
	
	for k, v in seqc.items():
		if v != max_model: return False
	
	return True


def entities_with_structures(entry):
	structloop = entry.get_loops_by_category("Atom_site")
	if len(structloop) != 1: return None
	else:                    structloop = structloop[0]
	
	model_id  = structloop.get_tag('Model_ID')
	entity_id = structloop.get_tag('Label_entity_ID')
	resid     = structloop.get_tag('Label_comp_ID')
	seqindex  = structloop.get_tag('Label_comp_index_ID')
	
	if not conformers_check(structloop): return None
	
	entseq = dict()
	prev = None
	counter = None
	for m, e, r, ind in zip(model_id, entity_id, resid, seqindex):
		if prev == None: 
			prev = m
			
		if m != prev: 
			for k,v in entseq.items():
				v = ''.join(v)
				entseq[k] = v
			return entseq
		
		if e not in entseq: entseq[e] = []
		
		if counter == None:
			if r not in aa_dict: return None
			entseq[e].append(aa_dict[r])
			counter = ind
			continue
		
		if ind != counter:
			if r not in aa_dict: return None
			entseq[e].append(aa_dict[r])
			counter = ind


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
		
		if seq_check(seq): return seq
		else: 			   return None
	
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
				if seq_check(seq): return seq
				else: 			   return None
		
		return None
		
	else:
		print(polymer[0])
		return None


def get_atom_types(aas, atms):
	
	atoms = list()
	
	for i in range(len(atms)):
		a = atms[i]
		if aas[i] not in aa_dict: continue
		if a not in atoms: 
			atoms.append(a)
	
	atoms.sort()
	return atoms


def get_shifts(shifts, proteins, id):
	eids = shifts.get_tag('Entity_ID')
	ind = shifts.get_tag('Seq_ID')
	aas = shifts.get_tag('Comp_ID')
	atm = shifts.get_tag('Atom_ID')
	val  = shifts.get_tag('Val')
	
	assert(len(eids) == len(ind) == len(aas) == len(atm) == len(val))
	
	shift_data = dict()
	atom_names = get_atom_types(aas,atm)
	
	for id in proteins.keys():
		shift_data[id] = dict()
		for an in atom_names:
			shift_data[id][an] = []
			for i in range(len(proteins[id])):
				shift_data[id][an].append(None)
	
	for i in range(len(eids)):
		id = eids[i]
		pos = int(ind[i])-1
		if aas[i] not in aa_dict:
			continue
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
		
		shift_data[id][atom][pos] = cs
	
	return shift_data


def get_coordinates(entry):
	# get save frame
	conf_set = entry.get_loops_by_category("Atom_site")
	
	model = conf_set.get_tag('Model_ID')
	entid = conf_set.get_tag('Label_entity_ID')
	seqid = conf_set.get_tag('Label_comp_index_ID')
	resid = conf_set.get_tag('Label_comp_ID')
	atmid = conf_set.get_tag('Label_atom_ID')
	cartx = conf_set.get_tag('Cartn_x')
	carty = conf_set.get_tag('Cartn_y')
	cartz = conf_set.get_tag('Cartn_z')
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
