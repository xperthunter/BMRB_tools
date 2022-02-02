#!/usr/bin/python3

import os
import sys
import re
import argparse
import json
from multiprocessing import Pool

import pynmrstar
import bmrb_tools as bt

parser = argparse.ArgumentParser(description='BMRB STAR3.0 to json')
parser.add_argument('--folder', required=True, type=str,
	metavar='<str>', help='File path to bmrb entries')
parser.add_argument('--json', required=True, type=str,
	metavar='<str>', help='Name of json file to dump to')
parser.add_argument('--err', required=True, type=str,
	metavar='<str>', help='file name for error output')
parser.add_argument('--csv', required=False, type=str,
	metavar='<str>', help='input csv file with specific ids to make json from')                           
parser.add_argument('--workers', required=False, type=int,
	metavar='<int>', help='number of workers for multiprocessing'
# grabbing all the args                                                          
arg = parser.parse_args()

def ent_parser(file):
	if not re.search(r'str$', file): return
	
	path = f'{arg.folder}{file}'
	if os.path.islink(path):
		path = os.path.join(os.path.dirname(arg.folder), os.readlink(path))
	
	try:
		ent = pynmrstar.Entry.from_file(path)
	except:
		return {'STAR': file}
	
	id = ent.get_tag('_Entry.ID')
	bid = id[0]
	
	if bt.assembly_number(ent) != 1:
		return {'assembly>1': file}
	
	shift_sets = bt.get_shiftlist(ent)
	
	if len(shift_sets) == 0:
		return {'!shiftloops': file}

	if len(shift_sets) != 1:
		return {'shiftlist!1': file}

	entities_with_cs = bt.get_ents_with_shifts(shift_sets[0])
	if len(entities_with_cs) is None:
		return {'ents_with_cs_error':file}
	
	seqs_per_ent = bt.seq_by_entity(ent,entities_with_cs)
	if seqs_per_ent is None:
		return {'seqsError':file}
	
	cs_data = bt.get_shifts(shift_sets[0], seqs_per_ent, bid)
	
	if cs_data is None:
		return {'shifterrors':file}
	
	for id in cs_data:
		dic = dict()
		dic = {
			'name' : bid,
			'seq' : seqs_per_ent[id]
		}
		
		for atm in cs_data[id].keys():
			dic[atm] = cs_data[id][atm]
	
	dic['structures'] = bt.get_coordinates(ent, dic['seq'])
	
	
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
	id_set = bt.csv_to_list(arg.csv)

data = list()

for file in os.listdir(arg.folder):
	if not re.search(r'str$', file): continue
	
	path = f'{arg.folder}{file}'
	if os.path.islink(path):
		path = os.path.join(os.path.dirname(arg.folder), os.readlink(path))
	
	try:
		ent = pynmrstar.Entry.from_file(path)
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
	
	if bt.assembly_number(ent) != 1:
		weird['assembly>1'].append(file)
		continue
	
	shift_sets = bt.get_shiftlist(ent)

	if len(shift_sets) == 0:
		weird['!shiftloops'].append(file)
		continue
	
	if len(shift_sets) != 1:
		weird['shiftlist!1'].append(file)
		continue
	
	entities_with_cs = bt.get_ents_with_shifts(shift_sets[0])
	if len(entities_with_cs) is None:
		weird['ents_with_cs_error'].append(file)
		continue
	
	seqs_per_ent = bt.seq_by_entity(ent,entities_with_cs)
	if seqs_per_ent is None:
		weird['seqsError'].append(file)
		continue
	
	cs_data = bt.get_shifts(shift_sets[0], seqs_per_ent, bid)
	
	if cs_data is None:
		weird['shifterrors'].append(file)
		continue
	
	for id in cs_data:
		dic = dict()
		dic = {
			'name' : bid,
			'seq' : seqs_per_ent[id]
		}
		
		for an in cs_data[id].keys():
			dic[an] = cs_data[id][an]
		
		data.append(dic)
	
	print('good')
	
with open(arg.json, 'w') as fp:
	json.dump(data, fp)
fp.close()

with open(arg.err, 'w') as fe:
	json.dump(weird, fe)
fe.close()