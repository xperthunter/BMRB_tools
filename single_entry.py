#!/usr/bin/python3

import os
import sys
import re
import argparse
import json

import pynmrstar
import bmrb_tools as bt


def entity_loops(ent):
	loops = ent.get_loops_by_category("Entity_assembly")
	if len(loops) > 1: return True
	else:              return False

def conf_loops(ent):
	loops = ent.get_loops_by_category("Atom_site")
	if len(loops) > 1: return True
	else:              return False

def shift_loops(ent):
	loops = ent.get_loops_by_category("Atom_chem_shift")
	if len(loops) > 1: return True
	else:              return False

def get_shiftframes(ent):
	loops = ent.get_saveframes_by_category("assigned_chemical_shifts")
	for loop in loops:
		yield loop

def conditions(condition_frame):
	type  = condition_frame.get_tag('Type')
	val   = condition_frame.get_tag('Val')
	units = condition_frame.get_tag('Val_units')
	
	for t,v,u in zip(type, val, units):
		yield t,v,u

"""
def entity_seqs(ent):
	entity_loop = ent.get_loops_by_category("Entity_assembly")[0]
	entids = entity_loop.get_tag("Entity_ID")
"""

file = sys.argv[1]

try:
	ent = pynmrstar.Entry.from_file(file)
	entry = dict()
	looped = dict()
except:
	print(f'could not process {file}')
	sys.exit()

id = ent.get_tag('_Entry.ID')
bid = id[0]
print(bid)
entry['name'] = bid

shifts_by_conditions = dict()

for frame in get_shiftframes(ent):
	frame_name      = frame['Sf_framecode'][0]
	condition_label = frame['Sample_condition_list_label'][0]
	
	condition_label = re.sub('^\$','',condition_label)
	if condition_label not in shifts_by_conditions: 
		shifts_by_conditions[condition_label] = list()
	
	shifts_by_conditions[condition_label].append(frame_name)

for k,v in shifts_by_conditions.items():
	print(k, v)




"""
for file in os.listdir(sys.argv[1]):
	print(file)
	#sys.exit()
	try:
		ent = pynmrstar.Entry.from_file(f'{sys.argv[1]}{file}')
	except:
		print(f'could not process {file}')
		continue
	
	id = ent.get_tag('_Entry.ID')
	bid = id[0]
	print(bid)
	
	if entity_loops(ent): continue
	if conf_loops(ent):   continue
	
	entity_loop = ent.get_loops_by_category("Entity_assembly")
	entids = entity_loop[0].get_tag("Entity_ID")
	print(entids)
	print(len(entids))
	if len(entids) > 1:
		if shift_loops(ent): 
			sloops = ent.get_loops_by_category("Atom_chem_shift")
			print(entity_loop[0])
			print(entids)
			print(sloops)
			print(file)
			sys.exit()
"""
"""
for file in os.listdir(sys.argv[1]):
	try:
		ent = pynmrstar.Entry.from_file(f'{sys.argv[1]}{file}')
	except:
		print(f'could not parse {file}')
#		sys.exit()

	shift_loops = tools.get_shiftlist(ent)
	if len(shift_loops) == 0: continue
	
	(num_comp,) = ent.get_tag('_Assembly.Number_of_components')
	
	if num_comp == '2':
		print(file)
#		sys.exit()
#	else:
#		print(file)
"""	
