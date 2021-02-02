"""
The original continuum_imaging_selfcal.py script is adapted for mosaics. It roughly 
follows the steps:
 1. Read input. Preamble.
 2. Split selfcal ms file. Same as continuum ms.
 3. Dirty image.
 4. Pre-selfcal image. 
 5. Populate model column in the ms file (tclean should do it in the previous step).
 6. Gaincal -> create gaintable list.
     6.5 Clearcal.
 7. Applycal gaintable list.
 8. Imaging. Increase iteration by 1. GOTO 5.
 9. Clearcal. Applycal list -> Final iteration clean.
10. Some comparisons pre-selfcal postselfcal, using the selcal model on the pre-selfcal... (?)

This script will more or less follow these.
"""
dryRun = True

# STAGE 1. PREAMBLE DEFINITIONS. Determination of array, antennae selection, et cetera.
# Import of several packages.
pre_preamble = set(globals())
import inspect
src_file_path = inspect.getfile(lambda: None)
execfile(os.path.join(os.path.dirname(src_file_path),"imaging_preamble.py"))

variables = ['imaging_parameters', 'from_cmd', 'image', 'almaimf_rootdir', 'exclude_7m', 
    'do_bsens', 'imaging_root', 'arrayname', 'logprint', 'opt_imaging_pars', 'metadata', 
    'bands', 'glob', 'imagingOptimalParameters', 'vis_image_parameters', 'fh', 'copy', 
    'continuum_files_per_field', 'only_7m', 'continuum_mses', 'check_output', 
    'dirtyImage', 'json_load_byteified', 'field_image_parameters']
variables.sort()
logprint("Defined variables in preamble: "+','.join(variables))

# STAGE 2. DEFINE AND/OR CREATE SELFCAL MS FILES.
from fAEG import create_selfcal_ms_files
selfcal_mses = create_selfcal_ms_files(mses=continuum_mses,vis_data=vis_image_parameters, arrayname=arrayname, do_bsens = do_bsens)



