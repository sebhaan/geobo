# Wrapper around main geobo script to include settings filename as argument.
# Copyright 2020 Sebastian Haan, The University of Sydney

import os 
import sys

if len(sys.argv) == 2:
	if os.path.isfile(sys.argv[1]):
		fname_settings = str(sys.argv[1])
	else:
		print('File ' + str(argv[1]) + ' does not exist')
		exit()
else:
	fname_settings ='settings.yaml'
	print("No settings file specified, using settings in: " + fname_settings)

from geobo import run_geobo