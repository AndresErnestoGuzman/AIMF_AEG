import inspect
src_file_path = inspect.getfile(lambda: None)
execfile(os.path.join(os.path.dirname(src_file_path),"imaging_preamble.py"))

bands = {metadata['ebs'][uid]['band'] for uid in metadata['ebs'].keys()}
fields = set(sum([metadata[band].keys() for band in bands],[]))
if 'FIELD_ID' in os.environ:
    fields = fields & set(os.getenv('FIELD_ID').split())

for band in bands:
	for field in fields:
		mymd = metadata[band][field]
		allspws = set(sum(mymd['spws'],[]))
		for spw in allspws:
			l_vises = []
			uids = []
			antenna = []
			for path,vis,spws in zip(mymd['path'],mymd['vis'],mymd['spws']):
				if spw in spws:
					uid = vis.replace('.ms.split.cal','')
					uids.append(uid)
					l_vis = os.path.join(path,'_'.join([uid,field,band,'spw'+str(spw)])+'.split')
					l_vises.append(l_vis)
					for vv in vis_image_parameters:
						if uid == vis_image_parameters[vv]['ebname']:
							antenna.append(vis_image_parameters[vv][arrayname+'_antennae'])
			for robust in [0]:
				imagename = os.path.join(imaging_root, "{0}_{1}_{2}_robust{3}_spw{4}".format(field,band,arrayname,robust,spw))
				tclean(vis = l_vises,
						field = field,
						spw = '0',
						antenna = antenna,
						imagename = imagename,
						imsize = field_image_parameters[band][field]['imsize'],
						cell = field_image_parameters[band][field]['cellsize'],
						phasecenter = field_image_parameters[band][field]['phasecenter'][0],
						specmode = 'cube',
						outframe = 'LSRK',
						chanchunks = -1,
						pblimit = 0.1,
						weighting = 'briggs',
						robust = robust,
						niter = 0,
						usemask = 'pb',
						pbmask = 0.25,
					)