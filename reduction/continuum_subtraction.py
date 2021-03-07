import os, json, sys
execfile("AIMF_AEG/reduction/defineRootdir.py")
from taskinit import mstool

from parse_contdotdat import parse_contdotdat, contchannels_to_linechannels # noqa: E402
from metadata_tools import logprint, json_load_byteified

ms = mstool()

with open('metadata.json', 'r') as fh:
    metadata = json_load_byteified(fh)
bands = {metadata['ebs'][uid]['band'] for uid in metadata['ebs'].keys()}
fields = set(sum([metadata[band].keys() for band in bands],[]))
if 'FIELD_ID' in os.environ:
    fields = fields & set(os.getenv('FIELD_ID').split())

for band in bands:
	for field in fields:
		mymd = metadata[band][field]
		for path,vis,spws in zip(mymd['path'],mymd['vis'],mymd['spws']):
			uid = vis.replace('.ms.split.cal','')
			mous = path.split('/')[-2]
			cont_sel = parse_contdotdat(mymd['contfit.dat'][mous])
			for spw in [3]:#spws:
				lfile = os.path.join(path,'_'.join([uid,field,band,'spw'+str(spw)])+'.split')
				if os.path.exists(lfile):
					ms.open(lfile)
				else:
					raise Exception("{0} does not exits!".format(lfile))
				
				freqs = {}
				freqs[0] = ms.cvelfreqs(spwids=[0],outframe='LSRK')
				ms.close()
				linechans = contchannels_to_linechannels(cont_sel, freqs)
				lfile_csub = lfile+'.contsub'
				if os.path.exists(lfile_csub):
					logprint("File {0}. Overwriting.".format(lfile_csub))
				uvcontsub(vis = lfile,
						field = field,
						fitspw = linechans,
						fitorder = 1,
						excludechans = True)
			
