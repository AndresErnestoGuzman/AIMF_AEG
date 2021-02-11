import os, json, time, sys
import numpy as np

from taskinit import casalog
from taskinit import msmdtool
from taskinit import mstool, tbtool
from tasks import split, flagmanager, flagdata, rmtables, concat

execfile("AIMF_AEG/reduction/defineRootdir.py")

from parse_contdotdat import parse_contdotdat, contchannels_to_linechannels # noqa: E402
from metadata_tools import logprint, json_load_byteified


msmd = msmdtool()
ms = mstool()
tb = tbtool()

# band name : frequency range (GHz)
bands = {'B3': (80, 110), 'B4': (125,163), 'B5': (163, 211),
         'B6': (211, 275), 'B7': (275, 373)
        }

# Function definitions. Logprint and some to get "strings" 
# from the metadata.json file. https://stackoverflow.com/a/33571117/1091
# Select fields with metadata. Combine with those defined in the FIELD_ID
# environmental variable
if True:    
    
    logprint("ALMAIMF_ROOTDIR directory set to {0}".format(os.getenv('ALMAIMF_ROOTDIR')))

    with open('metadata.json', 'r') as fh:
        metadata = json_load_byteified(fh)

    with open('contdatfiles.json', 'r') as fh:
        contdat_files = json.load(fh)

    # extract the fields from the metadata
    all_fields = set([])
    for band in bands:
        all_fields = all_fields | set(str(x) for x in metadata[band])

    if not os.getenv('FIELD_ID') == None:
        all_fields = all_fields & set(os.getenv('FIELD_ID').split())

    logprint("All fields include: {0}".format(", ".join(all_fields)))

    fields = set([])
    for band in bands:
        for field in all_fields:
            if field not in metadata[band]:
                logprint("Skipping {0}:{1} because it has no associated metadata"
                         .format(band, field))
            else:
                logprint("Processing continuum for {0}:{1} because it has metadata"
                         .format(band, field))
                fields = fields | set([field])
                metadata[band][field]['visfiles'] = []
                for path, vis in zip(metadata[band][field]['path'],metadata[band][field]['vis']):
                    metadata[band][field]['visfiles'].append(os.path.join(path,vis))

    logprint("Fields with metadata included in the 'fields' variable are: {0}".format(", ".join(fields)))
    
#  Determine contfiles. Split Best-Sensitivity continuum.
if True: 
    contdatFiles_assignment = {}
    for band in bands:
        contdatFiles_assignment[band] = {}
        for field in fields:
            contdatFiles_assignment[band][field] = {}
            # Define local metadata variable for the field and band
            if field not in metadata[band]:
                continue
            mymd = metadata[band][field]
            logprint("Metadata for {0}:{1} is {2}".format(band, field, mymd))
            for path, vis, spws, muid in zip(mymd['path'], mymd['vis'], mymd['spws'], mymd['muid']):
                contfile = mymd['cont.dat'][muid] 
                if os.path.exists(contfile):
                    contdatFiles_assignment[band][field][os.path.join(path,vis)] = contfile
                else:
                    logprint("****** No cont.dat file found for {0} = {1}:{2}.  "
                             .format(path, band, field))
                    logprint("No cont.dat file: Skipping - this file will not be included in the merged continuum.")

    
    visibilities_to_average_continuum = set([])
    
    for band in bands:
        for field in fields:
            for visfile in contdatFiles_assignment[band][field]:
                visibilities_to_average_continuum = visibilities_to_average_continuum | set([visfile])

    logprint("Visibilities for which we will flag line channels and split are: {0}"
                     .format(visibilities_to_average_continuum))

    for visi in visibilities_to_average_continuum:
        (visi_path, visi_filename) = os.path.split(visi)
        eb = re.sub('.*(uid[^.]+).*','\\1',visi_filename)
        fields_to_average_for_continuum = set(metadata['ebs'][eb]['fields']) & fields
        field = next(iter(fields_to_average_for_continuum))
        for band in bands:
            if field in metadata[band]:
                idx = metadata[band][field]['visfiles'].index(visi)
                widths = metadata[band][field]['continuum_widths_in_channels'][idx]
                spws = metadata[band][field]['spws'][idx]
                vis = metadata[band][field]['vis'][idx]
                path = metadata[band][field]['path'][idx]
                contvis_bestsens = str(os.path.join(path, "continuum_"+vis+"_bsens.cont"))
                tb.open(visi)
                if 'CORRECTED_DATA' in tb.colnames():
                    datacolumn = 'corrected'
                else:
                    datacolumn = 'data'
                tb.close()
                print "split(vis = {0}, spw = {1} ,field = {2}, outputvis = {3}, width = {4}, datacolumn = {5})".format(visi,
                    ",".join(map(str, spws)),",".join(list(fields_to_average_for_continuum)),contvis_bestsens,widths,datacolumn)
                 
                not_redo_contbsens = False
                if os.path.exists(contvis_bestsens):
                    not_redo_contbsens = msmd.open(contvis_bestsens)
                    if not_redo_contbsens:
                        not_redo_contbsens = not_redo_contbsens and fields_to_average_for_continuum.issubset(set(msmd.fieldnames()))
                        if not_redo_contbsens:
                            logprint("\nNot re-doing the BSENS continuum. File in existence seem to cover the requisites.")
                            logprint("To redo the bsens file please remove {0} from the calibrated folder.\n".format(contvis_bestsens))
                    else:
                        logprint("Continuum BSENS {0} not accesible. Removing and re-doing".format(contvis_bestsens))
                        rmtables(contvis_bestsens)
                        os.system("rm -rf "+contvis_bestsens)
                        os.system("rm -rf "+contvis_bestsens+".flagversions")
                    msmd.close()

                if not_redo_contbsens:
                    logprint("NOT doing continuum BSENS.")
                else:
                    if os.path.exists(contvis_bestsens):
                        logprint("Re-doing and replacing BSENS continuum file.")
                        rmtables(contvis_bestsens)
                        os.system("rm -rf "+contvis_bestsens)
                        os.system("rm -rf "+contvis_bestsens+".flagversions")
                    assert split(vis = visi,
                                 spw = ",".join(map(str, spws)),
                                 field = ",".join(list(fields_to_average_for_continuum)),
                                 outputvis = contvis_bestsens,
                                 width = widths,
                                 datacolumn = datacolumn), "Bsens Split Failed "
                    flagdata(vis = contvis_bestsens, mode = 'manual', autocorr = True, flagbackup = False) 
            else:
                logprint("No metadata associated with field {0} in {1}. Skipping.".format(field,band))
                    
#  Purest continuum line flags.
if True:
    for visi in visibilities_to_average_continuum:
        
        flagmanager(vis = visi, mode='delete', versionname='before_line_flags')
        flagmanager(vis = visi, mode='save'  , versionname='before_line_flags') # 
        (visi_path, visi_filename) = os.path.split(visi)
        eb = re.sub('.*(uid[^.]+).*','\\1',visi_filename)
        fields_to_average_for_continuum = set(metadata['ebs'][eb]['fields']) & fields
        
        ms.open(visi)
        
        for field in fields_to_average_for_continuum:
            for band in bands:
                if field in metadata[band]:
                    idx = metadata[band][field]['visfiles'].index(visi)
                    widths = metadata[band][field]['continuum_widths_in_channels'][idx]
                    spws = metadata[band][field]['spws'][idx]
                    vis = metadata[band][field]['vis'][idx]
                    path = metadata[band][field]['path'][idx]
                    frqlimsLSRK = metadata[band][field]['freqs'][idx]

                    try:
                        contdatfile = contdatFiles_assignment[band][field][visi]
                    except:
                        print("Problems in the definition of cont.dat file dictionary. contdatFiles_assignment = {0}".
                            format(contdatFiles_assignment))
                    
                    logprint("Using file {0} to define continuum in {1}.".format(contdatfile,visi))
                    cont_channel_selection = parse_contdotdat(contfile)
                    freqs = {}
                    for spw  in spws:
                        freqs[spw] = ms.cvelfreqs(spwids=[spw], outframe='LSRK')

                    linechannels, linefracs = contchannels_to_linechannels(cont_channel_selection,
                                                                freqs,
                                                                return_fractions=True)
                    logprint("Line fractions for {0}, field={1} are: {0}".format(vis, field, linefracs))
                    logprint("--> Flagging lines.\n")
                    flagdata(vis = visi, 
                            mode = 'manual', 
                            field = field, 
                            spw = linechannels, 
                            flagbackup = False)
                else:
                        logprint("No metadata associated with field {0} in {1}. Skipping.".format(field,band))

        ms.close()
        flagmanager(vis = visi, mode = 'delete', versionname = 'line_channel_flags')
        flagmanager(vis = visi, mode = 'save', versionname = 'line_channel_flags')
                    
# Splitting Purest continuum. Unflag lines in original files.
if True:
    cont_mses = []
    for visi in visibilities_to_average_continuum:
        (visi_path, visi_filename) = os.path.split(visi)
        eb = re.sub('.*(uid[^.]+).*','\\1',visi_filename)
        fields_to_average_for_continuum = set(metadata['ebs'][eb]['fields']) & fields
        field = next(iter(fields_to_average_for_continuum))
        for band in bands:
            if field in metadata[band]:
                idx = metadata[band][field]['visfiles'].index(visi)
                widths = metadata[band][field]['continuum_widths_in_channels'][idx]
                spws = metadata[band][field]['spws'][idx]
                vis = metadata[band][field]['vis'][idx]
                path = metadata[band][field]['path'][idx]
                
                contvis = str(os.path.join(path, "continuum_"+vis+".cont"))
                cont_mses.append(contvis) 
                
                tb.open(visi)
                if 'CORRECTED_DATA' in tb.colnames():
                    datacolumn = 'corrected'
                else:
                    datacolumn = 'data'
                tb.close()
                print "split(vis = {0}, spw = {1} ,field = {2}, outputvis = {3}, width = {4}, datacolumn = {5})".format(visi,",".
                    join(map(str, spws)),",".join(list(fields_to_average_for_continuum)),contvis,widths,datacolumn)
                
                contvis_exists = os.path.exists(contvis)
                not_redo_contbsens = False
                if contvis_exists:
                    not_redo_contbsens = msmd.open(contvis)
                    if not_redo_contbsens:
                        not_redo_contbsens = not_redo_contbsens and fields_to_average_for_continuum.issubset(set(msmd.fieldnames()))
                        if not_redo_contbsens:
                            logprint("\nNot re-doing the PUREST continuum. File in existence seem to cover the requisites.")
                            logprint("To redo the bsens file please remove {0} from the calibrated folder.\n".format(contvis))
                    else:
                        logprint("Continuum PUREST {0} not accesible. Removing and re-doing".format(contvis))
                        rmtables(contvis)
                        os.system("rm -rf "+contvis)
                        os.system("rm -rf "+contvis+".flagversions")
                    msmd.close()

                if not_redo_contbsens:
                    logprint("NOT doing continuum PUREST.")
                else:
                    if contvis_exists:
                        logprint("Re-doing and replacing PUREST continuum file.")
                    rmtables(contvis)
                    os.system("rm -rf "+contvis)
                    os.system("rm -rf "+contvis+".flagversions")
                    assert split(vis = visi,
                                 spw = ",".join(map(str, spws)),
                                 field = ",".join(list(fields_to_average_for_continuum)),
                                 outputvis = contvis,
                                 width = widths,
                                 datacolumn = datacolumn), "Purest Split Failed."
                    flagdata(vis = contvis, mode = 'manual', autocorr = True, flagbackup = False) 
            else:
                logprint("No metadata associated with field {0} in {1}. Skipping.".format(field,band))

    for visi in visibilities_to_average_continuum:
        flagmanager(vis = visi, mode='restore', versionname='before_line_flags')

with open('continuum_mses.txt', 'w') as fh:
    for line in cont_mses:
        fh.write(line+'\n')


with open('metadata_updated.json', 'w') as fh:
    json.dump(metadata, fh)

logprint("Completed split_cont_windows")
