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

from fAEG import create_selfcal_ms_files, selfcal_name, check_model_column, sc_iter_image
from imaging_parameters import selfcal_imaging_pars,  selfcal_pars

# STAGE 2. DEFINE AND/OR CREATE SELFCAL MS FILES.

selfcal_mses = create_selfcal_ms_files(mses=continuum_mses,vis_data=vis_image_parameters, eb_data=metadata['ebs'], arrayname=arrayname, do_bsens = do_bsens)
selfcal_files_per_field = {}
for band in bands:
    selfcal_files_per_field[band] = {}
    for field in continuum_files_per_field[band].keys():
            selfcal_files_per_field[band][field] = [selfcal_name(ms=cvis, arrayname=arrayname, field=field, band=band, bsens=False) for cvis in continuum_files_per_field[band][field]]

# STAGE 3. DIRTY IMAGE. Skip, use continuum_imaging_AEG.py

# STAGE 4 and 5. PRE-SELFCAL OR ITER 0. SHALLOW CLEAN.
sc_iter_image(iteracion=0,arrayname=arrayname,selfcal_files_per_field=selfcal_files_per_field,do_bsens=do_bsens,
    field_image_parameters=field_image_parameters,vis_image_parameters=vis_image_parameters,continuum_files_per_field=continuum_files_per_field,
    selfcal_imaging_pars=selfcal_imaging_pars,imaging_root=imaging_root,dryRun = dryRun, savemodel='modelcolumn',datacolumn='data')

#STAGE 5
caltables = {}
for sms in selfcal_mses:
    (directory,msfilename) = os.path.split(sms)
    caltables[msfilename] = []

for scp in selfcal_pars:
    match = re.match(r"([^_]+)_(B[1-9]+)_([127]+M)_.*robust([^\._]+)",scp)
    if match:
        (field, band, array, robust) = match.group(1,2,3,4)
    else:
        logprint("Unrecognized key in selfcal_pars", origin="cis_script")
        raise Exception()
    
    if not os.getenv('FIELD_ID') == None:
        if not field in set(os.getenv('FIELD_ID').split()):
            continue

    for iteracion in sorted(selfcal_pars[scp].keys()):
        if iteracion==0:
            continue
        # Next for is for all mses associated with the specific field
        for sms in selfcal_files_per_field[band][field]:
            (directory,msfilename) = os.path.split(sms)
            caltable = re.sub('.*(uid.*)_selfcal\.ms','\\1_iter'+str(iteracion)+'_'+str(selfcal_pars[scp][iteracion]['solint'])+'.cal',msfilename)
            if dryRun:
                logprint("Dry run. gaincal(vis={vis}, caltable={caltable},gaintable={gaintable}, {restofpars})".
                    format(vis = sms, caltable = caltable, gaintable = caltables[msfilename], restofpars=selfcal_pars[scp][iteracion]),origin="cis_script")
            else:
                logprint("gaincal. Iteration={0}".format(iteracion),origin="cis_script")
                gaincal(vis = sms,
                    caltable = caltable,
                    gaintable = caltables[msfilename],
                    **selfcal_pars[scp][iteracion])
            
            caltables[msfilename].append(caltable)

            if dryRun:
                logprint("Dry run. applycal(vis={vis}, gaintable={gaintable}),interp = 'linear', applymode = 'calonly', calwt = False)".
                    format(vis = sms,  gaintable = caltables[msfilename], origin = "cis_script"))
            else:
                logprint("applycal. Iteration={0}".format(iteracion), origin = "cis_script")
                clearcal(vis = sms, addmodel = True)
                applycal(vis = sms,
                        #gainfield = 
                        gaintable = caltables[msfilename],
                        interp = "linear",
                        applymode = 'calonly',
                        calwt = False)

        sc_iter_image(iteracion = iteracion, arrayname = arrayname, bands = set([band]), fields=set([field]), selfcal_files_per_field = selfcal_files_per_field, do_bsens = do_bsens, 
            field_image_parameters = field_image_parameters, vis_image_parameters = vis_image_parameters, continuum_files_per_field = continuum_files_per_field, 
            selfcal_imaging_pars = selfcal_imaging_pars, imaging_root = imaging_root, dryRun = dryRun, savemodel = 'modelcolumn', datacolumn = 'corrected')
           
