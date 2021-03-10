
dryRun = False
onlyDirty = False
continuum_subtracted = False
parallel = True
select_spws = {3}
robusts = [0]

fullwidth = '500km/s' # '450km/s'
lineFrequencies = {	'H41a':'92034.43415MHz','H39a':'106737.35656MHz','H29a':'256302.03519MHz',
					'H51b':'93607.31579MHz','H50b':'99225.20843MHz','H49b':'105301.85742MHz'}
# for 2016.1.00732.S 
# 'H41a' -> SPW1
# 'H39a' -> SPW3
line = 'H39a'
restfreq = lineFrequencies[line] if line else None

import inspect, glob
src_file_path = inspect.getfile(lambda: None)
execfile(os.path.join(os.path.dirname(src_file_path),"imaging_preamble.py"))

from fAEG import UniversalSet, rms_from_mad, line_imaging_dict, create_clean_model
if not select_spws:
	select_spws = UniversalSet()

from taskinit import msmdtool, qatool
msmd = msmdtool()
qa = qatool()

bands = {metadata['ebs'][uid]['band'] for uid in metadata['ebs'].keys()}
fields = set(sum([metadata[band].keys() for band in bands],[]))
if 'FIELD_ID' in os.environ:
    fields = fields & set(os.getenv('FIELD_ID').split())

suffix_cs = '.contsub' if continuum_subtracted else ''
suffix_cs_ima = '' if continuum_subtracted else '_noCsub'

for band in bands:
	for field in fields:
		mymd = metadata[band][field]
		allspws = set(sum(mymd['spws'],[])) & select_spws
		for spw in allspws:
			l_vises = []
			uids = []
			antenna = []
			for path,vis,spws in zip(mymd['path'],mymd['vis'],mymd['spws']):
				if spw in spws:
					uid = vis.replace('.ms.split.cal','')
					uids.append(uid)
					l_vis = os.path.join(path,'_'.join([uid,field,band,'spw'+str(spw)])+'.split'+suffix_cs)
					l_vises.append(l_vis)
					for vv in vis_image_parameters:
						if uid == vis_image_parameters[vv]['ebname']:
							antenna.append(vis_image_parameters[vv][arrayname+'_antennae'])
			for robust in robusts:
				imagename = os.path.join(imaging_root, "{0}_{1}_{2}_robust{3}_spw{4}{5}{6}".format(field,band,arrayname,robust,spw,suffix_cs_ima,line))
				
				msmd.open(l_vises[0]) # Assume in l_vises all have compatible channel widths in l_vises
				chw = np.mean(msmd.chanwidths(0))
				mean_freq_spw = qa.quantity(msmd.meanfreq(0),"Hz")
				msmd.close()
				if restfreq:
					chw_q = qa.convertdop(qa.div(qa.convert("{0}Hz".format(chw),'Hz'),qa.convert(qa.quantity(restfreq),'Hz')),'km/s')
				else:
					chw_q = qa.convertdop(qa.div(qa.convert("{0}Hz".format(chw),'Hz'),mean_freq_spw),'km/s')

				line_imaging_parameters = line_imaging_dict(chanwidth = qa.tos(chw_q),restfreq=restfreq,fullwidth=fullwidth)
				sys.exit()
				del line_imaging_parameters['width']
				

				if os.path.exists(imagename+".image"):
					logprint("{0}.image exists. Skipping dirty imaging.".format(imagename))
				elif dryRun:
					logprint("""tclean(vis = {vis}, field = {field},spw = {spw},antenna = {antenna},imagename = {imagename},
						imsize = {imsize} ,cell = {cell},phasecenter = {phasecenter}, specmode = {specmode},
						outframe = {outframe},chanchunks = {chanchunks},pblimit = {pblimit}, weighting={weighting},robust = {robust},niter = {niter},
						usemask = {usemask},pbmask = {pbmask},{rest})""".format(vis = l_vises,field = field,spw = '0',antenna = antenna,imagename = imagename,
							imsize = field_image_parameters[band][field]['imsize'],cell = field_image_parameters[band][field]['cellsize'],
							phasecenter = field_image_parameters[band][field]['phasecenter'][0],specmode = 'cube',
							outframe = 'LSRK',chanchunks = -1,pblimit = 0.1,weighting = 'briggs',robust = robust,
							niter = 0,usemask = 'pb',pbmask = 0.25,rest=line_imaging_parameters),origin='image_lines_script')
				else:
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
						pbmask = 0.25, parallel=parallel,
						**line_imaging_parameters
					)
					veloType = True if line else False
					if parallel:
						for tipo in ["image","residual","mask"]:
							ima = "{0}.{1}".format(imagename,tipo)
							if os.path.exists(ima):
								exportfits(imagename = ima,fitsimage=ima+".fits",overwrite=True, velocity=veloType)

if onlyDirty:
	logprint("Only dirty cubes.",origin='image_lines_script')
else:
	automasking_config_per_field = {}
	for field in fields:
		if False:
			automasking_config_per_field[field] = {'sidelobethreshold':2.0,'noisethreshold':5.00,'minbeamfrac':0.3,'lownoisethreshold':1.5,'negativethreshold':7.0}
		else: # Short
			automasking_config_per_field[field] = {'sidelobethreshold':2.0,'noisethreshold':4.25,'minbeamfrac':0.3,'lownoisethreshold':1.5,'negativethreshold':15.0}

	
	for band in bands:
		for field in fields:
			mymd = metadata[band][field]
			allspws = set(sum(mymd['spws'],[])) & select_spws
			for spw in allspws:
				l_vises = []
				uids = []
				antenna = []
				for path,vis,spws in zip(mymd['path'],mymd['vis'],mymd['spws']):
					if spw in spws:
						uid = vis.replace('.ms.split.cal','')
						uids.append(uid)
						l_vis = os.path.join(path,'_'.join([uid,field,band,'spw'+str(spw)])+'.split'+suffix_cs)
						l_vises.append(l_vis)
						for vv in vis_image_parameters:
							if uid == vis_image_parameters[vv]['ebname']:
								antenna.append(vis_image_parameters[vv][arrayname+'_antennae'])
				for robust in robusts:

					image_filename = "{0}_{1}_{2}_robust{3}_spw{4}{5}{6}".format(field,band,arrayname,robust,spw,suffix_cs_ima,line)
					imagename = os.path.join(imaging_root, image_filename)

					try:
						os.system("cp -r "+imagename+".image "+imagename+"_dirty.image")
					except:
						Exception("No dirty image.")

					msmd.open(l_vises[0])
					chw = np.mean(msmd.chanwidths(0))
					mean_freq_spw = qa.quantity(msmd.meanfreq(0),"Hz")
					msmd.close()
					if restfreq:
						chw_q = qa.convertdop(qa.div(qa.convert("{0}Hz".format(chw),'Hz'),qa.convert(qa.quantity(restfreq),'Hz')),'km/s')
					else:
						chw_q = qa.convertdop(qa.div(qa.convert("{0}Hz".format(chw),'Hz'),mean_freq_spw),'km/s')

					line_imaging_parameters = line_imaging_dict(chanwidth = qa.tos(chw_q),restfreq=restfreq,fullwidth=fullwidth)
					del line_imaging_parameters['width']
					
					if dryRun:
						logprint("""tclean(vis = {vis}, field = {field},spw = {spw},antenna = {antenna},imagename = {imagename},
							imsize = {imsize} ,cell = {cell},phasecenter = {phasecenter}, specmode = {specmode},
							outframe = {outframe},chanchunks = {chanchunks},pblimit = {pblimit}, weighting={weighting},robust = {robust},niter = {niter},
							usemask = {usemask},{aumc},{rest})""".format(
								vis = l_vises,field = field,spw = '0',antenna = antenna,imagename = imagename,
								imsize = field_image_parameters[band][field]['imsize'],cell = field_image_parameters[band][field]['cellsize'],
								phasecenter = field_image_parameters[band][field]['phasecenter'][0],specmode = 'cube',
								outframe = 'LSRK',chanchunks = -1,pblimit = 0.1,weighting = 'briggs',robust = robust,
								niter = int(1e5), usemask = 'auto-multithresh',aumc=automasking_config_per_field[field],
								rest=line_imaging_parameters),origin='image_lines_script')
					else:						
						if not continuum_subtracted:
							contimagename = sorted(glob.glob(os.path.join(imaging_root,"{field}_{band}_{arrayname}_robust{robust}_selfcal*.image.tt0".
								format(field=field,band=band,arrayname=arrayname,robust=robust))))[-1]
							
							if os.path.exists(contimagename):
								startmodel = create_clean_model(cubeimagename = imagename, contimagename = contimagename)
								if os.path.exists(imagename + ".model"):
									rmtables(imagename + ".model")
									os.system("rm -rf "+imagename + ".model")

								if parallel:
									wd = os.path.join(imaging_root,image_filename + ".workdirectory")
									if os.path.exists(wd):
										rmtables(os.path.join(wd,image_filename +".n*.model"))
										os.system("rm -rf "+os.path.join(wd,image_filename +"n*.model"))

								line_imaging_parameters['startmodel'] = startmodel

						line_imaging_parameters.update(automasking_config_per_field[field])	
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
							#chanchunks = -1,
							pblimit = 0.1,
							weighting = 'briggs',
							robust = robust,
							niter = int(1e5),
							nsigma = 3.0,
							usemask = 'auto-multithresh', pbmask=0.25,
							parallel = parallel,
							**line_imaging_parameters
						)
						if parallel:
							veloType = True if line else False
							for tipo in ["image","model","residual","mask"]:
								ima = "{0}.{1}".format(imagename,tipo)
								if os.path.exists(ima):
									exportfits(imagename = ima,fitsimage=ima+".fits",overwrite=True, velocity=veloType)

