import subprocess, os, glob

import os, json, sys
execfile("AIMF_AEG/reduction/defineRootdir.py")
from metadata_tools import logprint, json_load_byteified

if os.path.exists('metadata.json'):
	with open('metadata.json', 'r') as fh:
	    metadata = json_load_byteified(fh)
else:
	logprint("No metadata.json file. Ending execution", origin='cal_lines_script')

bands = {metadata['ebs'][uid]['band'] for uid in metadata['ebs'].keys()}
fields = set(sum([metadata[band].keys() for band in bands],[]))
if 'FIELD_ID' in os.environ:
    fields = fields & set(os.getenv('FIELD_ID').split())

#cube_vis = subprocess.check_output(['find','.','-name','*spw*split']).split()
ebs = set([])
l_vises = set([])

for band in bands:
	for field in fields:
		mymd = metadata[band][field]
		for path,vis,spws in zip(mymd['path'],mymd['vis'],mymd['spws']):
			uid = vis.replace('.ms.split.cal','')
			ebs.add(uid)
			caltables = []
			allcaltables = []
			allcaltables = glob.glob('_'.join([uid,field,band])+"*.cal")
			for cc in allcaltables:
				if not 'bsens' in cc:
					caltables.append(cc)
			caltables.sort()
			for spw in [0]:#spws:
				l_vis = os.path.join(path,'_'.join([uid,field,band,'spw'+str(spw)])+'.split')
				clearcal(vis = l_vis)
				applycal(vis = l_vis,
						gaintable = caltables,
				        interp = 'linear',
				        applymode = 'calonly',
				        calwt = False)



