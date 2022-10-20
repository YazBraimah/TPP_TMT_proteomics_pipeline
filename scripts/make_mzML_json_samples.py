#!/usr/bin/env python3
'''
Make a samples.json file with sample names and file names.
'''
def msg(name=None):                                                            
    return ''' make_json_samples.py <samples_files>

    mzML file names should have the following format:
        <sample_name>.mzML
        '''

import json
from glob import glob
from sys import argv
import sys
import argparse
import glob
parser = argparse.ArgumentParser(description='Make a samples.json file with sample names and file names.', usage=msg())


mzmlFiles = argv[1:]
mzMLs = []
for mz in mzmlFiles:
    mzMLs.extend(glob.glob(mz))

FILES = {}

# Change this line to extract a sample name from each filename.
SAMPLES = [mzML.split('/')[-1].split('.')[0] for mzML in mzMLs]

for sample in SAMPLES:
    file1 = lambda mzML: sample in mzML and 'mzML' in mzML
    if any('mzXML' in s for s in sorted(filter(file1, mzMLs))):
       
        FILES[sample] = {}
        FILES[sample]= sorted(filter(file1, mzMLs))
    else:
        
        FILES[sample] = {}
        FILES[sample]= sorted(filter(file1, mzMLs))



js = json.dumps(FILES, indent = 4, sort_keys=True)
open('samples_mzML.json', 'w').writelines(js)
