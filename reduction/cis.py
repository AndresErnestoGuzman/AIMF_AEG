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

from fAEG import create_selfcal_ms_files, selfcal_name, check_model_column, sc_iter_image, rms_from_mad
from imaging_parameters import selfcal_imaging_pars,  selfcal_pars

# STAGE 2. DEFINE AND/OR CREATE SELFCAL MS FILES.

selfcal_mses = create_selfcal_ms_files(mses=continuum_mses,vis_data=vis_image_parameters, 
    eb_data=metadata['ebs'], arrayname=arrayname, do_bsens = do_bsens)
selfcal_files_per_field = {}
for band in bands:
    selfcal_files_per_field[band] = {}
    for field in continuum_files_per_field[band].keys():
            selfcal_files_per_field[band][field] = [selfcal_name(ms=cvis, arrayname=arrayname, field=field, band=band, bsens=False) for cvis in continuum_files_per_field[band][field]]

# STAGE 3. DIRTY IMAGE. Skip, use continuum_imaging_AEG.py

# STAGE 4 and 5. PRE-SELFCAL OR ITER 0. SHALLOW CLEAN. DELETE MODEL COLUMN OF VISIBILITIES
for dbs in do_bsens:
    sc_iter_image(iteracion=0,arrayname=arrayname,selfcal_files_per_field=selfcal_files_per_field,do_bsens=dbs,
        field_image_parameters=field_image_parameters,vis_image_parameters=vis_image_parameters,continuum_files_per_field=continuum_files_per_field,
        selfcal_imaging_pars=selfcal_imaging_pars,imaging_root=imaging_root,dryRun = dryRun, savemodel='modelcolumn',datacolumn='data')

#STAGE 5
caltables = {}
for sms in selfcal_mses:
    (directory,msfilename) = os.path.split(sms)
    caltables[msfilename] = []

def filtrar(st):
    match = re.match(r"([^_]+)_(B[1-9]+)_([127]+M)_.*robust([^\._]+)",st)
    return "_".join(match.group(1,2,3,4))
auxdict = {v1:v2 for v1,v2 in zip(map(filtrar,selfcal_pars.keys()),selfcal_pars.keys())}
filtered_keys = auxdict.keys()

for scp in filtered_keys:
    match = re.match(r"([^_]+)_(B[1-9]+)_([127]+M)_([^\._]+)",scp)
    if match:
        (field, band, array, robust) = match.group(1,2,3,4)
        scp = re.sub(r'_([^_]+)$','_robust\\1',scp)
    else:
        raise Exception("Unrecognized key in selfcal_pars")
    
    if not os.getenv('FIELD_ID') == None:
        if not field in set(os.getenv('FIELD_ID').split()):
            continue

    for dbs in do_bsens:
        for iteracion in sorted(selfcal_pars[scp].keys()): # just iterations. Same number bsens and purest.
            if iteracion==0:
                continue
            # Next for loop is for all mses associated with the specific field
            for sms in selfcal_files_per_field[band][field]:            
                
                key = scp
                suffix = '_bsens' if dbs else ''
                sms = sms.replace("selfcal.ms","selfcal"+suffix+".ms")
                (directory,msfilename) = os.path.split(sms)     
                
                caltable = re.sub('.*(uid.*)_selfcal.*','\\1_iter'+str(iteracion)+str(selfcal_pars[key][iteracion]['calmode'])+
                    '_'+str(selfcal_pars[key][iteracion]['solint'])+suffix+'.cal',msfilename)

                if os.path.exists(caltable):
                    logprint("Iteration={0}. Calibration table {1} already exist!. Using that one. Skipping this calibration iteration.".
                        format(iteracion,caltable),origin="cis_script")
                    continue

                if dryRun:
                    logprint("Dry run. gaincal(vis={vis}, caltable={caltable},gaintable={gaintable}, {restofpars})".
                        format(vis = sms, caltable = caltable, gaintable = caltables[msfilename], restofpars=selfcal_pars[key][iteracion]),origin="cis_script")
                else:
                    logprint("gaincal. Iteration={0}".format(iteracion),origin="cis_script")
                    gaincal(vis = sms,
                        caltable = caltable,
                        gaintable = caltables[msfilename],
                        **selfcal_pars[key][iteracion])
                
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

            imagenames = sc_iter_image(iteracion = iteracion, arrayname = arrayname, bands = set([band]), fields=set([field]), selfcal_files_per_field = selfcal_files_per_field,
                do_bsens = dbs, field_image_parameters = field_image_parameters, vis_image_parameters = vis_image_parameters, continuum_files_per_field = continuum_files_per_field, 
                selfcal_imaging_pars = selfcal_imaging_pars, imaging_root = imaging_root, dryRun = dryRun, savemodel = 'modelcolumn', datacolumn = 'corrected')

            if len(imagenames)>1:
                raise Exception("imagenames should be size one list in this case")

            
            # MOAR CLEAN in LAST ITERATION
            if iteracion == max(map(int,sorted(selfcal_pars[scp].keys()))) and not dryRun:
                logprint("MOAR CLEAN.", origin='cis_script')
                last_iter_imaging_pars = {scp: selfcal_imaging_pars[scp]}
                rms = rms_from_mad(imagenames[0]+".image.tt0")
                last_iter_imaging_pars[scp]['threshold'][iteracion] = "{0:.2f}mJy".format(2500*rms)
                sc_iter_image(iteracion = iteracion, arrayname = arrayname, bands = set([band]), fields=set([field]), selfcal_files_per_field = selfcal_files_per_field,
                    do_bsens = dbs, field_image_parameters = field_image_parameters, vis_image_parameters = vis_image_parameters, 
                    continuum_files_per_field = continuum_files_per_field, selfcal_imaging_pars = last_iter_imaging_pars, imaging_root = imaging_root, 
                    dryRun = dryRun, savemodel = 'none', datacolumn = 'corrected')

            logprint("Done imaging of {0}, iteration={1} of {2}".format(imagenames[0],iteration, max(map(int,sorted(selfcal_pars[scp].keys())))), origin="cis_scrip")