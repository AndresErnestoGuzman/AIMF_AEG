#imaging

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

import inspect
src_file_path = inspect.getfile(lambda: None)
# https://stackoverflow.com/questions/16771894/python-nameerror-global-name-file-is-not-defined
execfile(os.path.join(os.path.dirname(src_file_path),"defineRootdir.py"))

#from getversion import git_date, git_version
from metadata_tools import logprint, json_load_byteified
#from make_custom_mask import make_custom_mask

from imaging_parameters import imaging_parameters
from fAEG import imagingOptimalParameters, dirtyImage, image

#from tasks import tclean, exportfits, plotms, split
#from taskinit import  iatool
#from utils import validate_mask_path

if os.path.exists('metadata.json'):
    with open('metadata.json', 'r') as fh:
        metadata = json_load_byteified(fh)
else:
    raise FileError("metadata.json not found!")


#msmd = msmdtool()
#ia = iatool()

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
if os.getenv('DO_BSENS').lower() == 'true':
    do_bsens.append(True)
    logprint("Using BSENS measurement set")
if os.getenv('DO_BSENS_ONLY').lower() == 'true':
    do_bsens = [True]

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
