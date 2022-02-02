#!/usr/bin/python3

import argparse
import json 
import os
import re
import sys

import pynmrstar
import bmrb_tools as bt

parser = argparse.ArgumentParser(description='BMRB STAR3.0 to json')
parser.add_argument('--folder', required=True, type=str,
	metavar='<str>', help='File path to bmrb entries')
parser.add_argument('--json', required=True, type=str,
	metavar='<str>', help='Name of json file to dump to')
parser.add_argument('--err', required=True, type=str,
	metavar='<str>', help='file name for error output')
#parser.add_argument('--workers')

arg = parser.parse_args()

weird = {
	'STAR' : [],
	'assembly>1' : [],
	'!shiftloops' : [],
	'shiftlist!1' : [],
	'seqsError' : [],
	'shifterrors' : [],
	'ents_with_cs_error' : [],
	'not_requested' : [],
	'no_conditions': [],
	'conformer_issue': []
}
data = list()
for file in os.listdir(arg.folder):
	# only STAR files in folder
	if not re.search(r'str$', file): continue
	
	path = f'{arg.folder}{file}'
	# check if softlink
	if os.path.islink(path):
		path = os.path.join(os.path.dirname(arg.folder), os.readlink(path))
	
	# try and read .str with pynmrstar
	try:
		ent = pynmrstar.Entry.from_file(path)
	except:
		weird['STAR'].append(file)
		continue
	print(file)
	
	id = ent.get_tag('_Entry.ID')
	bid = id[0]
	
	# dont want str files with more than 1 assembly
	if bt.assembly_number(ent) != 1:
		weird['assembly>1'].append(file)
		continue
	
	shift_sets = bt.get_shiftlist(ent)
	
	# must have a chemical shift loop
	if len(shift_sets) == 0:
		weird['!shiftloops'].append(file)
		continue
	
	# if more than 1 shift loop, skip 
	if len(shift_sets) != 1:
		weird['shiftlist!1'].append(file)
		continue
	
	# get sample conditions
	sample_conditions = bt.conditions(ent)
	if sample_conditions == None: # if getting conditions did not work
		weird['no_conditions'].append(file)
		continue
	
	print(json.dumps(sample_conditions,indent=2))
	
	ews = bt.entities_with_structures(ent)
	if ews == None: 
		weird['confomer_issue'].append(file)
	print(json.dumps(ews,indent=2))
	
	sys.exit()
	
		
"""
identify sample conditions X
entities in conformer models X
entity assemly -> entities 
	what entities are in the assembly
	for each obtain their sequence
	entities in assembly not in structure
	entities with shifts
	entities with shifts and structures have same sequence
atom types with shifts
empty shift lists
fill shifts lists and shift dictionaries
conformers xyz paried with shifts
"""