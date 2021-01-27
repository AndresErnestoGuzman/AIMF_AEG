import os, json, time, sys
import numpy as np

from taskinit import casalog
from taskinit import msmdtool
from taskinit import mstool, tbtool
from tasks import split, flagmanager, flagdata, rmtables, concat

if 'almaimf_rootdir' in locals():
    os.environ['ALMAIMF_ROOTDIR'] = almaimf_rootdir

if os.getenv('ALMAIMF_ROOTDIR') is None:
    try:
        import metadata_tools
        os.environ['ALMAIMF_ROOTDIR'] = os.path.split(metadata_tools.__file__)[0]
    except ImportError:
        raise ValueError("metadata_tools not found on path; make sure to "
                         "specify ALMAIMF_ROOTDIR environment variable "
                         "or your PYTHONPATH variable to include the directory"
                         " containi ng the ALMAIMF code.")
else:
    sys.path.append(os.getenv('ALMAIMF_ROOTDIR'))


msmd = msmdtool()
ms = mstool()
tb = tbtool()

# band name : frequency range (GHz)
bands = {'B3': (80, 110),
         'B6': (210, 250),
        }

# Function definitions. Logprint and some to get "strings" 
# from the metadata.json file. https://stackoverflow.com/a/33571117/10916251
if True:
    def logprint(string):
        casalog.post(string, origin='split_cont_windows')
        print(string)

    def ebMatch(absfilename,uid,path):
        return (absfilename in uid) and (os.path.split(absfilename)[0] == path)

    def json_load_byteified(file_handle):
        return _byteify(
            json.load(file_handle, object_hook=_byteify),
            ignore_dicts=True
        )

    def json_loads_byteified(json_text):
        return _byteify(
            json.loads(json_text, object_hook=_byteify),
            ignore_dicts=True
        )

    def _byteify(data, ignore_dicts = False):
        # if this is a unicode string, return its string representation
        if isinstance(data, unicode):
            return data.encode('utf-8')
        # if this is a list of values, return list of byteified values
        if isinstance(data, list):
            return [ _byteify(item, ignore_dicts=True) for item in data ]
        # if this is a dictionary, return dictionary of byteified keys and values
        # but only if we haven't already byteified it
        if isinstance(data, dict) and not ignore_dicts:
            return {
                _byteify(key, ignore_dicts=True): _byteify(value, ignore_dicts=True)
                for key, value in data.iteritems()
            }
        # if it's anything else, return it in its original form
        return data


# Select fields with metadata. Combine with those defined in the FIELD_ID
# environmental variable
if True:    
    
    logprint("ALMAIMF_ROOTDIR directory set to {0}".format(os.getenv('ALMAIMF_ROOTDIR')))

    with open('metadata.json', 'r') as fh:
        metadata = json_load_byteified(fh)

    with open('contdatfiles.json', 'r') as fh:
        contdat_files = json.load(fh)

    # extract the fields from the metadata
    all_fields = set(str(x) for x in metadata['B3']) | set(str(x) for x in metadata['B6'])
    if not os.getenv('FIELD_ID') == None:
        all_fields = all_fields & set(os.getenv('FIELD_ID').split())

    logprint("all_fields include: {0}".format(", ".join(all_fields)))

    fields = set([])
    for band in bands:
        for field in all_fields:
            if field not in metadata[band]:
                logprint("Skipping {0}:{1} because it has no metadata"
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
                logprint("Skipping {0}:{1} because it has no metadata"
                         .format(band, field))
                continue
            else:
                logprint("Processing continuum for {0}:{1} because it has metadata"
                         .format(band, field))

            mymd = metadata[band][field]
            logprint("Metadata for {0}:{1} is {2}".format(band, field, mymd))
            for path, vis, spws, muid in zip(mymd['path'], mymd['vis'], mymd['spws'], mymd['muid']):
                contfile = mymd['cont.dat'][muid] 
                if not os.path.exists(contfile):
                    logprint("****** No cont.dat file found for {0} = {1}:{2}.  "
                             .format(path, band, field))
                    if field + band + muid in contdat_files:
                        # raise ValueError("Going to use a different cont.dat for this config?")
                        contfile = contdat_files[field + band + muid]
                        if not os.path.exists(contfile):
                            logprint("No cont.dat file: Skipping - this file will not be included in the merged continuum.")
                            continue
                        else:
                            logprint("No cont.dat file: Using {0} instead.".format(contfile))
                    else:
                        logprint("No cont.dat file: Skipping - this file will not be included in the merged continuum.")
                        continue
                
                visibilityFile = os.path.join(path,vis)
                contdatFiles_assignment[band][field][visibilityFile] = contfile

    
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
                print "split(vis = {0}, spw = {1} ,field = {2}, outputvis = {3}, width = {4}, datacolumn = {5})".format(visi,",".
                    join(map(str, spws)),",".join(list(fields_to_average_for_continuum)),contvis_bestsens,widths,datacolumn)
                
                contvis_bestsens_exists = os.path.exists(contvis_bestsens)
                not_redo_contbsens = False
                if contvis_bestsens_exists:
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
                    if contvis_bestsens_exists:
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
                    for spw, spwlim in zip(spws,frqlimsLSRK):
                        freqs[spw] = spwlim
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


# if False:
#     for band in bands:

#         for field in fields:
#         mymd = metadata[band][field]
        
#         for path, vis, spws, muid, frqlimsLSRK, widths in zip(mymd['path'], mymd['vis'], mymd['spws'], mymd['muid'], mymd['freqs'],mymd['continuum_widths_in_channels']):
#             visibilityFile = os.path.join(path,vis)
#             if visibilityFile in contdatFiles_assignment[band][field]:
#                 freqs = {}
#                 for spw, spwlim in zip(spws,frqlimsLSRK);
#                     freqs[spw] = spwlim

#                 contvis = str(os.path.join(path, "continuum_"+vis+".cont"))
#                 contvis_bestsens = str(os.path.join(path, "continuum_"+vis+"_bsens.cont"))


#                 contfile = contdatFiles_assignment[band][field][visibility]
#                 cont_channel_selection = parse_contdotdat(contfile)
#                 linechannels, linefracs = contchannels_to_linechannels(cont_channel_selection,
#                                                                 freqs,
#                                                                 return_fractions=True)
#                 logprint("Line fractions for {0}, field={1} are: {0}".format(vis, field, linefracs))
#                 # not clear why this is done in other imaging scripts, but it
#                 # seems to achieve the wrong effect.
#                 # initweights(vis=visfile, wtmode='weight', dowtsp=True)
#                 flagdata(vis=visibilityFile, mode='manual', field = field, spw=linechannels, flagbackup=False)
#                 flagmanager(vis=visibilityFile, mode='save', versionname='line_channel_flags')


#     visibilities_postlineflag = visibilities_to_save_preflag_state
#     for visi in visibilities_postlineflag:
#         flagmanager(vis=str(visi), mode='save', versionname='line_channel_flags')

#     for band in bands:

#         for field in fields:
#         mymd = metadata[band][field]
#         for path, vis, spws, frqlimsLSRK, widths in zip(mymd['path'], mymd['vis'], mymd['spws'], mymd['freqs'],mymd['continuum_widths_in_channels']):
#             visibilityFile = os.path.join(path,vis)
#             if visibilityFile in contdatFiles_assignment[band][field]:


#     #sys.exit()

#


# sys.exit()

# cont_mses = []
# cont_mses_unconcat = []

# # split the continuum data
# for band in bands:
#     cont_to_merge[band] = {}
#     for field in all_fields:

#         cont_to_merge[band][field] = []

#         # Check whether there is metadata for this field and band. 
#         # If not, skip to next field
#         if field not in metadata[band]:
#             logprint("Skipping {0}:{1} because it has no metadata"
#                      .format(band, field))
#             continue
#         else:
#             logprint("Processing continuum for {0}:{1} because it has metadata"
#                      .format(band, field))

#         # Define local metadata variable for the field and band
#         mymd = metadata[band][field]
#         logprint("Metadata for {0}:{1} is {2}".format(band, field, mymd))

#         # Iteration through the metadata entries. Goes through all visibilities
#         # which include this field and band.
#         for path, vis, spws, muid in zip(mymd['path'], mymd['vis'], mymd['spws'], mymd['muid']):

#             t0 = time.time() # Measuring the time it takes to run this script
            
#             # The specific cont.dat file for this visibility. Test for existence either as a file or as part of the
#             # special dictionary contdat_files
#             contfile = mymd['cont.dat'][muid] 
#             if not os.path.exists(contfile):
#                 logprint("****** No cont.dat file found for {0} = {1}:{2}.  "
#                          .format(path, band, field))
#                 if field + band + muid in contdat_files:
#                     # raise ValueError("Going to use a different cont.dat for this config?")
#                     contfile = contdat_files[field + band + muid]
#                     if not os.path.exists(contfile):
#                         logprint("No cont.dat file: Skipping - this file will not be included in the merged continuum.")
#                         continue
#                     else:
#                         logprint("No cont.dat file: Using {0} instead.".format(contfile))
#                 else:
#                     logprint("No cont.dat file: Skipping - this file will not be included in the merged continuum.")
#                     continue
            

#             # Channel selection of the continuum. Commiting to use contfile and
#             # put information in the metadata dictionaries
#             cont_channel_selection = parse_contdotdat(contfile)
#             mymd['cont.dat_file'] = contfile
#             contdat_files[field + band + muid] = contfile

#             # Define continuum visibility. 
#             visfile = str(os.path.join(path, vis))
#             contvis = str(os.path.join(path, "continuum_"+vis+".cont"))
#             contvis_bestsens = str(os.path.join(path, "continuum_"+vis+"_bsens.cont"))

#             cont_to_merge[band][field].append(contvis)

#             if os.path.exists(contvis) and os.path.exists(contvis_bestsens):
#                 logprint("Skipping width determination for {0} = {1}:{2} because "
#                          "it's done (both for bsens & cont)".format(contvis, band, field),)
#             else:
#                 logprint("Determining widths for {0} to {1}, {2}:{3}"
#                          .format(visfile, contvis, band, field),)

#                 # determine target widths
#                 msmd.open(visfile)
#                 ms.open(visfile)
#                 Synth_HPBW = 0.3 #Smallest synth HPBW among target sample in arcsec
#                 # values interpolated by Roberto from https://science.nrao.edu/facilities/vla/docs/manuals/oss2016A/performance/fov/bw-smearing
#                 PB_HPBW = 21. * (300. / bands[band][0]) # PB HPBW at lowest band freq
#                 # targetwidth = 10e6 # 10 MHz
#                 targetwidth = 0.25 * (Synth_HPBW / PB_HPBW) * bands[band][0] * 1e9 # 98% BW smearing criterion
#                 widths = []
#                 freqs = {}
#                 logprint("Determining smoothing widths for continuum data.  "
#                          "PB_HPBW = {0}, targetwidth = {1}".format(PB_HPBW, targetwidth))
#                 for spw in spws:
#                     chwid = np.abs(np.mean(msmd.chanwidths(spw)))
#                     wid = int(targetwidth/chwid)
#                     if wid <= 0:
#                         logprint("Failure at chwid = {0}, wid = {1}.  ".format(chwid, wid))
#                         raise ValueError("The channel width is greater than "
#                                          "the target line width for spw {0} "
#                                          "in ms {1}".format(spw, visfile))
#                     else:# 
#                         # CASA *cannot* handle wid > nchan
#                         # This one also insists that there will be at least 3
#                         # output channels in all cases. Also, divides nchan.
#                         wid = int(msmd.nchan(spw) / 3)
#                         divisors = np.array([1])
#                         _nn = msmd.nchan(spw)
#                         for _pp in [2,3,5]:
#                             _prodaux = np.outer(divisors,_pp**np.where(_nn%_pp**np.arange(0,int(ceil(log(_nn,_pp))))==0)[0])
#                             divisors = np.sort(np.ndarray.flatten(_prodaux))
#                         divisors = divisors[np.where((divisors < wid) & (divisors > _nn/3))]
#                         if len(divisors)>0:
#                             wid = divisor.max()

#                     widths.append(wid)
#                     # these are TOPO freqs: freqs[spw] = msmd.chanfreqs(spw)
#                     try:
#                         freqs[spw] = ms.cvelfreqs(spwid=[spw], outframe='LSRK')
#                     except TypeError:
#                         freqs[spw] = ms.cvelfreqs(spwids=[spw], outframe='LSRK')

#                 msmd.close()
#                 ms.close()

#             if not os.path.exists(contvis) or not os.path.exists(contvis_bestsens):
#                 tb.open(visfile)
#                 if 'CORRECTED_DATA' in tb.colnames():
#                     datacolumn = 'corrected'
#                 else:
#                     datacolumn = 'data'
#                 tb.close()


#             if os.path.exists(contvis):
#                 logprint("Continuum: Skipping {0} because it's done".format(contvis),)
#             elif field not in fields:
#                 logprint("Skipping {0} because it is not one of the "
#                          "selected fields (but its metadata is being "
#                          "collected in continuum_mses.txt)".format(contvis))
#             else:
#                 logprint("Flagging and splitting {0} to {1} for continuum"
#                          .format(visfile, contvis),)
#                 logprint("contfile is {0}".format(contfile))


#                 linechannels, linefracs = contchannels_to_linechannels(cont_channel_selection,
#                                                             freqs,
#                                                             return_fractions=True)
#                 logprint("Line fractions are: {0}".format(linefracs))


#                 flagmanager(vis=visfile, mode='save',
#                             versionname='before_cont_flags')

#                 # not clear why this is done in other imaging scripts, but it
#                 # seems to achieve the wrong effect.
#                 # initweights(vis=visfile, wtmode='weight', dowtsp=True)


#                 flagdata(vis=visfile, mode='manual', spw=linechannels,
#                          flagbackup=False)


#                 flagmanager(vis=visfile, mode='save',
#                             versionname='line_channel_flags')

#                 rmtables(contvis) # redundant. contvis does not exist
#                 os.system('rm -rf ' + contvis + '.flagversions')


#                 # Average the channels within spws
#                 # (assert here checks that this completes successfully)
#                 assert split(vis=visfile,
#                              spw=",".join(map(str, spws)),
#                              field=field,
#                              outputvis=contvis,
#                              width=widths,
#                              datacolumn=datacolumn), "Split failed!"

#                 if not os.path.exists(contvis):
#                     raise IOError("Split failed for {0}".format(contvis))

#                 # If you flagged any line channels, restore the previous flags
#                 flagmanager(vis=visfile, mode='restore',
#                             versionname='before_cont_flags')

#                 # flag out the autocorres
#                 flagdata(vis=contvis, mode='manual', autocorr=True)



#             if os.path.exists(contvis_bestsens):
#                 logprint("Skipping {0} because it's done".format(contvis_bestsens),)
#             elif field not in fields:
#                 logprint("Skipping {0} because it is not one of the "
#                          "selected fields (but its metadata is being "
#                          "collected in continuum_mses.txt)".format(contvis_bestsens))
#             else:
#                 logprint("Splitting 'best-sensitivity' {0} to {1} for continuum"
#                          .format(visfile, contvis_bestsens),)

#                 # Average the channels within spws for the "best sensitivity"
#                 # continuum, in which nothing is flagged out
#                 assert split(vis=visfile,
#                              spw=",".join(map(str, spws)),
#                              field=field,
#                              outputvis=contvis_bestsens,
#                              width=widths,
#                              datacolumn=datacolumn), "Split Failed 2"
#                 flagdata(vis=contvis_bestsens, mode='manual', autocorr=True)

#             logprint("Finished splitting for {0} to {1}, {2}:{3} in {4} seconds"
#                      .format(visfile, contvis, band, field, time.time() - t0))


#         member_uid = path.split("member.")[-1].split("/")[0]
#         merged_continuum_fn = os.path.join(path,
#                                            "{field}_{band}_{muid}_continuum_merged.cal.ms"
#                                            .format(field = field,
#                                                    band = band,
#                                                    muid = member_uid)
#                                           )

#         # merge the continuum measurement sets to ease bookkeeping
#         if os.path.exists(merged_continuum_fn):
#             logprint("Skipping merged continuum {0} because it's done"
#                      .format(merged_continuum_fn),)
#         elif field not in fields:
#             logprint("Skipping {0} because it is not one of the "
#                      "selected fields (but its metadata is being "
#                      "collected in continuum_mses.txt)".format(merged_continuum_fn))
#         else:
#             logprint("Merging continuum for {0} {1} into {2}"
#                      .format(merged_continuum_fn, field, band),)
#             if(len(cont_to_merge[band][field])>1):
#                 concat(vis = cont_to_merge[band][field],concatvis = merged_continuum_fn)
#                 flagdata(vis=merged_continuum_fn, mode='manual', autocorr=True)
#             else:
#                 os.system("mv "+cont_to_merge[band][field][0]+" "+merged_continuum_fn)

#         cont_mses.append(merged_continuum_fn)


#         # merge the best sensitivity continuum too
#         merged_continuum_bsens_fn = os.path.join(
#             path,
#             "{field}_{band}_{muid}_continuum_merged_bsens.cal.ms"
#             .format(field=field, band=band, muid=member_uid)
#         )

#         if os.path.exists(merged_continuum_bsens_fn):
#             logprint("Skipping merged continuum bsens {0} because it's done"
#                      .format(merged_continuum_bsens_fn),)
#         elif field not in fields:
#             logprint("Skipping {0} because it is not one of the "
#                      "selected fields (but its metadata is being "
#                      "collected in continuum_mses.txt)".format(merged_continuum_bsens_fn))
#         else:
#             logprint("Merging bsens continuum for {0} {1} into {2}"
#                      .format(merged_continuum_bsens_fn, field, band),)

#             # Note this search-and-replace pattern: we use this instead
#             # of separately storing the continuum bsens MS names
#             if(len(cont_to_merge[band][field])>1):
#                 concat(vis=[x.replace(".cont", "_bsens.cont")
#                             for x in cont_to_merge[band][field]],
#                        concatvis = merged_continuum_bsens_fn,)
#                 flagdata(vis=merged_continuum_fn, mode='manual', autocorr=True)
#             else:
#                 os.system("mv "+cont_to_merge[band][field][0].replace(".cont", "_bsens.cont")+
#                             " "+merged_continuum_bsens_fn)

#         # for debug purposes, we also track the split, unmerged MSes
#         cont_mses_unconcat += cont_to_merge[band][field]


with open('continuum_mses.txt', 'w') as fh:
    for line in cont_mses:
        fh.write(line+'\n')

# with open('continuum_mses_unconcat.txt', 'w') as fh:
#     for line in cont_mses_unconcat:
#         fh.write(line+'\n')

# with open('cont_metadata.json', 'w') as fh:
#     json.dump(cont_to_merge, fh)

with open('metadata_updated.json', 'w') as fh:
    json.dump(metadata, fh)

# with open('contdatfiles_updated.json', 'w') as fh:
#     json.dump(contdat_files, fh)

logprint("Completed split_cont_windows")
