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
        lines are flagged out.
    FIELD_ID=<name>
        If this parameter is set, filter out the imaging targets and only split
        fields with this name (e.g., "W43-MM1", "W51-E", etc.).
        Metadata will still be collected for *all* available MSes.
    BAND_TO_IMAGE=B3 or B6
        If this parameter is set, only image the selected band.

The environmental variable ``ALMAIMF_ROOTDIR`` should be set to the directory
containing this file.
"""

onlyDirtyImaging = False

import os, sys, argparse, re, glob, copy

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
        from metadata_tools import determine_imsize, determine_phasecenter, logprint, json_load_byteified
        os.environ['ALMAIMF_ROOTDIR'] = os.path.split(metadata_tools.__file__)[0]
    except ImportError:
        raise ValueError("metadata_tools not found on path; make sure to "
                         "specify ALMAIMF_ROOTDIR environment variable "
                         "or your PYTHONPATH variable to include the directory"
                         " containing the ALMAIMF code.")
else:
    import sys
    sys.path.append(os.getenv('ALMAIMF_ROOTDIR'))

almaimf_rootdir = os.getenv('ALMAIMF_ROOTDIR')



from getversion import git_date, git_version
from metadata_tools import determine_imsize, determine_phasecenter, logprint, json_load_byteified
#from make_custom_mask import make_custom_mask

from imaging_parameters import imaging_parameters
# Default sigma calculation for threshold
images_for_sigma_estimation = glob.glob("imaging_results/*dirty.image.tt0")
for i in images_for_sigma_estimation:
    (path,filename) = os.path.split(i)
    auxiliar = filename.split('_')
    field = auxiliar.pop(0)
    band = auxiliar.pop(0)
    array = auxiliar.pop(0)
    robust_value = re.sub('robust([^\.]+).*','\\1',auxiliar.pop(0))
    key = "{0}_{1}_{2}_robust{3}".format(field, band, array, robust_value)
    if not 'threshold' in imaging_parameters[key]:
        stats = imstat(i)
        ss = 1.4826*stats['medabsdevmed'][0]
        imaging_parameters[key]['threshold'] = {0: "{0:.2f}mJy".format(1000*2*ss)}


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


# load the list of continuum MSes from a file
# (this file has one continuum MS full path, e.g. /path/to/file.ms, per line)
with open('continuum_mses.txt', 'r') as fh:
    continuum_mses = [x.strip() for x in fh.readlines()]

do_bsens = [False]

if os.getenv('DO_BSENS') is not None and os.getenv('DO_BSENS').lower() != 'false':
    do_bsens.append(True)
    logprint("Using BSENS measurement set")
    # continuum_mses += [x.replace('_continuum_merged.cal.ms',
    #                              '_continuum_merged_bsens.cal.ms')
    #                    for x in continuum_mses]

if exclude_7m:
    arrayname = '12M'
elif only_7m:
    arrayname = '7M'
else:
    arrayname = '7M12M'

""" 
Populate dictionaries 
	continuum_files_per_field : The keys are the fields. Each field associated with a list of continuum ms files
	vis_image_parameters : best image parameters per visibility file
"""
continuum_files_per_field = {}
vis_image_parameters = {}
bands = set()
for continuum_ms in continuum_mses:
    # get EB name
    ebname =  re.sub('.*(uid[^.]+).*','\\1',os.path.split(continuum_ms)[1])
    band = metadata['ebs'][ebname]['band']
    bands = bands | {band}
    vis_image_parameters[continuum_ms] = {}
    vis_image_parameters[continuum_ms]['ebname'] = ebname

    # allow optional cmdline args to skip one or the other band
    if os.getenv('BAND_TO_IMAGE'):
        if band not in os.getenv('BAND_TO_IMAGE'):
            logprint("Skipping band {0} because it is not in {1}"
                     .format(band, os.getenv('BAND_TO_IMAGE')),
                     origin='almaimf_cont_imaging')
            continue
        logprint("Imaging only band {0}".format(os.getenv('BAND_TO_IMAGE')),
                 origin='almaimf_cont_imaging')
    if not band in continuum_files_per_field:
        continuum_files_per_field[band] = {}
    fields = metadata['ebs'][ebname]['fields']
    for field in fields:
        if field in continuum_files_per_field[band]:
            continuum_files_per_field[band][field].append(continuum_ms)
        else:
            continuum_files_per_field[band][field] = [continuum_ms]

    field = fields[0]
    coosys,racen,deccen = determine_phasecenter(ms = continuum_ms, field = field)
    phasecenter = "{0} {1}deg {2}deg".format(coosys, racen, deccen)
    vis_image_parameters[continuum_ms]['phasecenter'] = phasecenter
    (dra,ddec,pixscale) = list(determine_imsize(ms=continuum_ms, field=field,
                                                phasecenter=(racen,deccen),
                                                exclude_7m=exclude_7m,
                                                only_7m=only_7m,
                                                spw='all',
                                                pixfraction_of_fwhm=1/8. if only_7m else 1/4.))
    imsize = [dra, ddec]
    vis_image_parameters[continuum_ms]['imsize'] = imsize
    vis_image_parameters[continuum_ms]['cellsize'] = ['{0:0.2f}arcsec'.format(pixscale)] * 2
    vis_image_parameters[continuum_ms]['cellsize_in_arcsec'] = [floor(pixscale*1000 + 0.5)/1000] * 2  # in arcsec
    msmd.open(continuum_ms)
    antennae = ",".join([x for x in msmd.antennanames() if 'CM' not in x])
    vis_image_parameters[continuum_ms]['12M_antennae'] = antennae
    antennae = ",".join([x for x in msmd.antennanames() if 'CM' in x])
    vis_image_parameters[continuum_ms]['7M_antennae'] = antennae
    vis_image_parameters[continuum_ms]['7M12M_antennae'] = ""
    msmd.close()

# Determine best image parameters per field. You need to check all image parameters of each continuum file containing the field.
field_image_parameters = {}
bands  = continuum_files_per_field.keys()
for band in bands:
    fields = continuum_files_per_field[band].keys()
    field_image_parameters[band] = {}
    for field in fields:
        field_image_parameters[band][field] = {}
        field_image_parameters[band][field]['imsize'] = [0,0]
        field_image_parameters[band][field]['cellsize_in_arcsec'] = [1e6,1e6]
        field_image_parameters[band][field]['phasecenter'] = []
        for cvis in continuum_files_per_field[band][field]:
            field_image_parameters[band][field]['imsize'] = max(field_image_parameters[band][field]['imsize'],vis_image_parameters[cvis]['imsize'])
            field_image_parameters[band][field]['cellsize_in_arcsec'] = min(field_image_parameters[band][field]['cellsize_in_arcsec'],
                                                                    vis_image_parameters[cvis]['cellsize_in_arcsec'])
            field_image_parameters[band][field]['cellsize'] = ['{0:0.2f}arcsec'.format(field_image_parameters[band][field]['cellsize_in_arcsec'][0])]*2
            field_image_parameters[band][field]['phasecenter'].append(vis_image_parameters[cvis]['phasecenter'])


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
                if (not os.path.exists(imname+".image.tt0")):
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
                if onlyDirtyImaging:
                    logprint("Ony dirty images this run.")
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

