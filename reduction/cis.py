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
dryRun = False

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

from fAEG import create_selfcal_ms_files, selfcal_name, check_model_column

# STAGE 2. DEFINE AND/OR CREATE SELFCAL MS FILES.

selfcal_mses = create_selfcal_ms_files(mses=continuum_mses,vis_data=vis_image_parameters, eb_data=metadata['ebs'], arrayname=arrayname, do_bsens = do_bsens)
selfcal_files_per_field = {}
for band in bands:
    selfcal_files_per_field[band] = {}
    for field in continuum_files_per_field[band].keys():
            selfcal_files_per_field[band][field] = [selfcal_name(ms=cvis, arrayname=arrayname, field=field, band=band, bsens=False) for cvis in continuum_files_per_field[band][field]]

# STAGE 3. DIRTY IMAGE. Skip, use continuum_imaging_AEG.py

# STAGE 4 and 5. PRE-SELFCAL
for band in bands:
    fields = set(selfcal_files_per_field[band].keys())
    if not os.getenv('FIELD_ID') == None:
        fields = fields & set(os.getenv('FIELD_ID').split())
    logprint("Doing images for fields {0}".format(", ".join(list(fields))), origin='cis_preselfcal')
    for field in fields:
        for dbs in do_bsens:
            suffix = "_bsens" if dbs else ""
            phasecenter = field_image_parameters[band][field]['phasecenter'][0]
            antennae = [vis_image_parameters[vis][arrayname+'_antennae'] for vis in continuum_files_per_field[band][field]]
            selfcal_visibilities = [x.replace(".ms",suffix+".ms") for x in selfcal_files_per_field[band][field]]
            imsize = field_image_parameters[band][field]['imsize']
            cell = field_image_parameters[band][field]['cellsize']
            #contimagename = os.path.join(imaging_root, field) +"_"+ band +"_" + arrayname + suffix

            for robust in [0]:
                imname = os.path.join(imaging_root, field) +"_"+ band +"_" + arrayname + suffix+"_robust{0}_preselfcal".format(robust)
                if os.path.exists(imname):
                    logprint("Image {0} already exists. Skip redoing pre-selfcal image.".format(os.path.split(imname)[1]),origin='cis_preselfcal')
                else:
                    image(visibilities=selfcal_visibilities, cell=cell, field=field, band=band, arrayname=arrayname, robust=robust, imaging_root=imaging_root,
                            imsize=imsize, antennae = antennae, phasecenter=phasecenter, suffix=suffix, imaging_parameters = imaging_parameters, dryrun = dryRun,
                            savemodel='modelcolumn',datacolumn='data',imagename=imname)

                if not dryRun:
                    if all(check_model_column(vis=selfcal_visibilities)):
                        logprint("Model column populated from pre-sefcal step.",origin='cis_check_model_column')




