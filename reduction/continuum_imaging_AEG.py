"""
Continuum imaging scripts.  There must be a ``continuum_mses.txt`` file in the
directory this is run in.  That file is produced by ``split_windows.py``.impars

You can set the following environmental variables for this script:
    EXCLUDE_7M=<boolean>
        If this parameter is set (to anything), the 7m data will not be
        included in the images if they are present.
    DO_BSENS=False
        By default, the "best sensitivity" continuum images will also be made
        if the _bsens.ms measurement sets exist.  If you set DO_BSENS=False,
        they will be skipped instead.  "Best sensitivity" means that no spectral
        lines are flmsmdagged out.
    FIELD_ID=<name>
        If this parameter is set, filter out the imaging targets and only split
        fields with this name (e.g., "W43-MM1", "W51-E", etc.).
        Metadata will still be collected for *all* available MSes.
    BAND_TO_IMAGE=B3 or B6
        If this parameter is set, only image the selected band.

The environmental variable ``ALMAIMF_ROOTDIR`` should be set to the directory
containing this file.


Additional Notes
================
USE_SELFCAL_MS is an environmental variable you can set if you want the imaging
to be done using the selfcal.ms file instead of the default continuum MS file.
It is primarily for debug purposes and you shouldn't need it.
"""

onlyDirtyImaging = False
dryRun = True

import os, sys, argparse, re, glob, copy 
from subprocess import check_output

try:
    # If run from command line
    aux = os.path.dirname(os.path.realpath(sys.argv[2]))
    if os.path.isdir(aux) and os.path.exists(os.path.join(aux, 'continuum_imaging.py')):
        almaimf_rootdir = aux
        from_cmd = True
    else:
        from_cmd = False
except:
    from_cmd = False

if 'almaimf_rootdir' in locals():
    os.environ['ALMAIMF_ROOTDIR'] = almaimf_rootdir
if os.getenv('ALMAIMF_ROOTDIR') is None:
    try:
        import metadata_tools
        from metadata_tools import determine_phasecenter, logprint, json_load_byteified
        os.environ['ALMAIMF_ROOTDIR'] = os.path.split(metadata_tools.__file__)[0]
    except ImportError:
        raise ValueError("metadata_tools not found on path; make sure to "
                         "specify ALMAIMF_ROOTDIR environment variable "
                         "or your PYTHONPATH variable to include the directory"
                         " containing the ALMAIMF code.")
elif not os.getenv('ALMAIMF_ROOTDIR') in sys.path:
    sys.path.append(os.getenv('ALMAIMF_ROOTDIR'))

almaimf_rootdir = os.getenv('ALMAIMF_ROOTDIR')


from getversion import git_date, git_version
from metadata_tools import determine_phasecenter, logprint, json_load_byteified
#from make_custom_mask import make_custom_mask

from imaging_parameters import imaging_parameters
from fAEG import imagingOptimalParameters


from tasks import tclean, exportfits, plotms, split
from taskinit import msmdtool, iatool
from utils import validate_mask_path

if os.path.exists('metadata.json'):
    with open('metadata.json', 'r') as fh:
        metadata = json_load_byteified(fh)
else:
    raise FileError("metadata.json not found!")


msmd = msmdtool()
ia = iatool()

imaging_root = "imaging_results"
if not os.path.exists(imaging_root):
    os.mkdir(imaging_root)

# Command line options
if from_cmd:
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', nargs=1, 
            help='Casa parameter')
    parser.add_argument('--exclude7M', action='store_true',
            help='Include 7M data')
    parser.add_argument('--only7M', action='store_true',
            help='Only image 7M data')
    parser.add_argument('--use_selfcal', action='store_true',
            help='Use self-calibrated ms')
    args = parser.parse_args()
    exclude_7m = args.exclude7M
    only_7m = args.only7M
    if args.use_selfcal:
        os.environ['USE_SELFCAL_MS'] = "true"

if 'exclude_7m' not in locals():
    if os.getenv('EXCLUDE_7M') is not None:
        exclude_7m = bool(os.getenv('EXCLUDE_7M').lower() == 'true')
    else:
        exclude_7m = False

if 'only_7m' not in locals():
    if os.getenv('ONLY_7M') is not None:
        only_7m = bool(os.getenv('ONLY_7M').lower() == 'true')
    else:
        only_7m = False

do_bsens = [False]
if os.getenv('DO_BSENS') is not None and os.getenv('DO_BSENS').lower() != 'false':
    do_bsens.append(True)
    logprint("Using BSENS measurement set")

if exclude_7m:
    arrayname = '12M'
elif only_7m:
    arrayname = '7M'
else:
    arrayname = '7M12M'

# load the list of continuum MSes from a file
# (this file has one continuum MS full path, e.g. /path/to/file.ms, per line)
with open('continuum_mses.txt', 'r') as fh:
    continuum_mses = [x.strip() for x in fh.readlines()]

opt_imaging_pars = imagingOptimalParameters(continuum_mses,metadata, exclude_7m  = exclude_7m, only_7m = only_7m)
vis_image_parameters = opt_imaging_pars['visBased']
field_image_parameters = opt_imaging_pars['fieldBased']
continuum_files_per_field = opt_imaging_pars['filesPerField']

bands  = continuum_files_per_field.keys()

# The imaging
for band in bands:
    fields = set(continuum_files_per_field[band].keys())
    if not os.getenv('FIELD_ID') == None:
        fields = fields & set(os.getenv('FIELD_ID').split())
    logprint("Doing images for fields {0}".format(", ".join(list(fields))))
    for field in fields:
        for dbs in do_bsens:
            suffix = "_bsens" if dbs else ""
            phasecenter = field_image_parameters[band][field]['phasecenter'][0]
            antennae = [vis_image_parameters[cvis][arrayname+'_antennae'] for cvis in continuum_files_per_field[band][field]]
            continuum_visibilities = [x.replace("split.cal","split.cal"+suffix) for x in continuum_files_per_field[band][field]]
            imsize = field_image_parameters[band][field]['imsize']
            cell = field_image_parameters[band][field]['cellsize']
            contimagename = os.path.join(imaging_root, field) +"_"+ band +"_" + arrayname + suffix

            for robust in [0]:
                key = "{0}_{1}_{2}_robust{3}".format(field, band,arrayname,robust)
                impars = imaging_parameters[key]
                #print impars
                #sys.exit()
                impars = copy.copy(impars)
                dirty_impars = copy.copy(impars)
                if 'maskname' in dirty_impars:
                    del dirty_impars['maskname']
                dirty_impars['niter'] = 0
                dirty_impars['usemask'] = 'pb'
                dirty_impars['pbmask'] = 0.25
                if 'threshold' in dirty_impars:
                    del dirty_impars['threshold']

                ## Dirty imaging
                imname = contimagename+"_robust{0}_dirty".format(robust)
                if (not os.path.exists(imname+".image.tt0") and not dryRun):
                    logprint("Dirty imaging file {0}".format(imname),
                             origin='almaimf_cont_imaging')
                    #antennae = [vis_image_parameters[cvis][arrayname+'_antennae'] for cvis in continuum_visibilities]
                    tclean(vis = continuum_visibilities,
                           field = field.encode(),
                           imagename = imname,
                           #phasecenter = phasecenter,
                           outframe = 'LSRK',
                           veltype = 'radio',
                           #usemask='pb',
                           interactive = False,
                           cell = cell,
                           imsize = imsize,
                           antenna = antennae,
                           pbcor = True,
                           **dirty_impars
                          )

                    stats = imstat(imagename=imname+".image.tt0",algorithm='hinged-fences',fence=2)
                    mad = stats['medabsdevmed'][0]
                    rms  = 1.4826* mad
                    imhead(imname+".image.tt0",mode='put',hdkey='RMS',hdvalue=rms)
    
                    ia.open(imname+".residual.tt0")
                    ia.sethistory(origin='almaimf_cont_imaging',
                                  history=["{0}: {1}".format(key, val) for key, val in
                                           impars.items()])
                    ia.sethistory(origin='almaimf_cont_imaging',
                                  history=["git_version: {0}".format(git_version),
                                           "git_date: {0}".format(git_date)])
                    ia.close()
                    exportfits(imname+".image.tt0", imname+".image.tt0.fits", overwrite = True)
                else:
                	logprint("Skipping dirty image block. Image={image}, visibilities={vis}".format(image=imname,vis=continuum_visibilities))
                
                if onlyDirtyImaging:
                    logprint("Only dirty images this run.")
                    continue
 
                if not 'maskname' in impars:
                    impars['usemask'] = 'pb'
                    maskname = ''
                    if not 'pblimit' in impars:
                        raise KeyError('Either define maskname or set pblimit > 0')

                # for compatibility w/self-calibration: if a list of parameters is used,
                # just use the 0'th iteration's parameters
                
                
                for iteracion in impars['threshold'].keys():
                    impars_thisiter = copy.copy(impars)
                    
                    for key, val in impars_thisiter.items():
                        if isinstance(val, dict):
                            impars_thisiter[key] = val[iteracion]

                    if impars_thisiter['usemask'] == 'pb'and 'maskname' in impars_thisiter:
                        del impars_thisiter['maskname']
                    elif impars_thisiter['usemask'] == 'user':
                        if 'maskname' in impars_thisiter:
                            maskname = validate_mask_path(impars_thisiter['maskname'],
                                                          os.getenv('ALMAIMF_ROOTDIR'))
                            del impars_thisiter['maskname']
                        if 'mask' not in impars_thisiter:
                            impars_thisiter['mask'] = maskname

                    imname = contimagename+"_robust{0}".format(robust)

                    if not os.path.exists(imname+".image.tt0"):
                        logprint("Cleaning file {0} for the first time. New image.".format(imname),
                                 origin='almaimf_cont_imaging')
                    else:
                        logprint("Cleaning already existing file {0}. Will start from this image in tclean.".format(imname),
                                 origin='almaimf_cont_imaging')
                    
                    antennae = [vis_image_parameters[cvis][arrayname+'_antennae'] for cvis in continuum_files_per_field[band][field]]
                    # print """ tclean(vis = {0},
                    #        field = {1},
                    #        imagename = {2},
                    #        phasecenter = {3},
                    #        outframe = 'LSRK',
                    #        veltype = 'radio',
                    #        interactive = False,
                    #        cell = {4},
                    #        imsize = {5},
                    #        antenna = {6},
                    #        pbcor = True,
                    #        **impars_thisiter
                    #       )""".format(continuum_visibilities,field.encode(),imname,phasecenter,cellsize,imsize,antennae)
                    # sys.exit()
                    rmtables(imname+".mask")
                    if not dryRun:
	                    if os.path.exists(imname+".mask"):
	                        os.system("rm -rf "+imname+".mask")
	                    tclean(vis = continuum_visibilities,
	                           field = field.encode(),
	                           imagename = imname,
	                           #phasecenter = phasecenter,
	                           outframe = 'LSRK',
	                           veltype = 'radio',
	                           #usemask='user',
	                           interactive = False,
	                           cell = cell,
	                           imsize = imsize,
	                           antenna = antennae,
	                           pbcor = True,
	                           **impars_thisiter
	                          )

	                    stats = imstat(imagename=imname+".image.tt0",algorithm='hinged-fences',fence=2)
	                    mad = stats['medabsdevmed'][0]
	                    rms  = 1.4826* mad
	                    imhead(imname+".image.tt0",mode='put',hdkey='RMS',hdvalue=rms)
	 
	                    #os.system("cp -r {0}".image.tt0" {0}.iter{1}".format(imname,iteracion))
	                    ia.open(imname+".image.tt0")
	                    ia.sethistory(origin='almaimf_cont_imaging',
	                                  history=["{0}: {1}".format(key, val) for key, val in
	                                           impars.items()])
	                    ia.sethistory(origin='almaimf_cont_imaging',
	                                  history=["git_version: {0}".format(git_version),
	                                           "git_date: {0}".format(git_date)])
	                    ia.close()

	                    exportfits(imname+".image.tt0", imname+".image.tt0.fits", overwrite = True)
	                    exportfits(imname+".image.tt0.pbcor", imname+".image.tt0.pbcor.fits", overwrite = True)
	                    exportfits(imname+".model.tt0", imname+".model.tt0.fits", overwrite = True)
	                    exportfits(imname+".residual.tt0", imname+".residual.tt0.fits", overwrite = True)
                    else:
                    	logprint("Dry run: skipping imaging of {image} from visibilities {vis}.".format(image=imname,vis=continuum_visibilities))
