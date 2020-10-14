#!/usr/bin/python3

import os
import sys
import re
import argparse
import json

import pynmrstar
import bmrb2json_tools as tools

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
	