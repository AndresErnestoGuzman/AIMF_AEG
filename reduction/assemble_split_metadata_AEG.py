import os, glob, json, time, re, sys, subprocess
import numpy as np

from taskinit import casalog
from taskinit import msmdtool
from taskinit import mstool, tbtool


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
                         " containing the ALMAIMF code.")
elif not os.getenv('ALMAIMF_ROOTDIR') in sys.path:
    sys.path.append(os.getenv('ALMAIMF_ROOTDIR'))

from parse_contdotdat import parse_contdotdat, contchannels_to_linechannels # noqa: E402

msmd = msmdtool()
ms = mstool()
tb = tbtool()

# band name : frequency range (GHz)
bands = {'B3': (80, 110), 'B4': (125,163), 'B5': (163, 211),
         'B6': (211, 275), 'B7': (275, 373)
        }
antenna_diameters = {'7M':7, '12Mshort':12, '12Mlong':12}

def logprint(string):
    casalog.post(string, origin='assemble_ms_metadata')
    print(string)

logprint("ALMAIMF_ROOTDIR directory set to {0}".format(os.getenv('ALMAIMF_ROOTDIR')))

metadata = {b: {} for b in bands}
metadata['ebs'] = {}
contdat_files = {}

science_goals = glob.glob("science_goal*")

t1 = time.time()

for sg in science_goals:
        
    splitcalfiles = subprocess.check_output(['find',sg,'-name','*split.cal']).split()
    
    for filename in splitcalfiles:

        dirpath = "/".join(filename.split("/")[0:-1])
        
        t0 = time.time()
        logprint("Collecting metadata for {0}; t0={1}".format(filename, t0))
        
        msmd.open(filename)
        ms.open(filename)
        mssummary = ms.summary()
        ms.close()
        
        antnames = msmd.antennanames()
        fieldnames = np.array(msmd.fieldnames())
        fields = np.unique(fieldnames[msmd.fieldsforintent('OBSERVE_TARGET#ON_SOURCE')])

        eb = re.sub('.*(uid[^.]+).*','\\1',os.path.split(filename)[1])
        metadata['ebs'][eb] = {}
        metadata['ebs'][eb]['path'] = os.path.split(filename)[0]
        metadata['ebs'][eb]['file'] = os.path.split(filename)[1]
        metadata['ebs'][eb]['fields'] = fields.tolist()
        if mssummary['timeref'] == 'UTC':
            metadata['ebs'][eb]['beginTimeString'] = au.mjdToDatestring(mssummary['BeginTime'])
        else:
            sys.exit()

        #assert len(np.unique(field)) == 1, "ERROR: field={0} fieldnames={1}".format(field, fieldnames)
        for field in fields:
        
            frq0 = msmd.chanfreqs(0)
            for bb, (lo, hi) in bands.items():
                #print lo, hi
                try:
                    if lo*1e9 < frq0 and hi*1e9 > frq0:
                        band = bb
                except ValueError:
                    if lo*1e9 < np.min(frq0) and hi*1e9 > np.max(frq0):
                        band = bb
            
            if not 'band' in metadata['ebs'][eb]:
                metadata['ebs'][eb]['band'] = band

            if any('PM' in nm for nm in antnames):
                if len(antnames) <= 4:
                    with open(os.path.join(filename, "{0}_{1}_TP".format(field, band)), 'w') as fh:
                        fh.write("{0}".format(antnames))
                        logprint("Skipping total power MS {0} [Elapsed: {1}]".format(filename, time.time()-t0))
                        msmd.close()
                        continue
                else:
                    logprint("WARNING: MS {0} contains PM antennae but is apparently not a TP data set".format(filename))

            try:
                summary = msmd.summary()
            except RuntimeError:
                logprint("Skipping FAILED MS {0} [Elapsed: {1}]".format(filename, time.time() - t0))
                msmd.close()
                continue

            logprint("NOT skipping ms {0} [Elapsed: {1}]".format(filename, time.time() - t0))
            spws = msmd.spwsforfield(field)
            targetspws = msmd.spwsforintent('OBSERVE_TARGET*')
            # this is how DOSPLIT in scriptForPI decides to split
            spws = [int(ss) for ss in spws if (ss in targetspws) and (msmd.nchan(ss) > 4)]

            # muid is 1 level above calibrated
            muid = filename.split("/")[-3]

            # need the full ms to get LSRK frequencies
            ms.open(filename)
            try:
                frqs = [ms.cvelfreqs(spwid=[spw], outframe='LSRK') for spw in spws]
                frqdict = {spw: ms.cvelfreqs(spwid=[spw], outframe='LSRK') for spw in spws}
            except TypeError:
                frqs = [ms.cvelfreqs(spwids=[spw], outframe='LSRK') for spw in spws]
                frqdict = {spw: ms.cvelfreqs(spwids=[spw], outframe='LSRK') for spw in spws}
            ms.close()

            frqslims = [(frq.min(), frq.max()) for frq in frqs]

            if field in metadata[band]:
                metadata[band][field]['path'].append(os.path.abspath(dirpath)),
                metadata[band][field]['vis'].append(filename.split("/")[-1])
                metadata[band][field]['spws'].append(spws)
                metadata[band][field]['freqs'].append(frqslims)
                metadata[band][field]['muid'].append(muid)
            else:
                metadata[band][field] = {'path': [os.path.abspath(dirpath)],
                                         'vis': [filename.split("/")[-1]],
                                         'spws': [spws],
                                         'freqs': [frqslims],
                                         'muid': [muid],
                }

            for ii in range(0,len(mssummary.keys())):
                fnum = 'field_'+str(ii)
                if fnum in mssummary.keys():
                    fname = mssummary[fnum]['name']
                    if fname == field:
                        rarad = mssummary[fnum]['direction']['m0']['value']
                        decrad = mssummary[fnum]['direction']['m1']['value']
                        refsystem = mssummary[fnum]['direction']['refer']
                        if not(refsystem=='J2000' or refsystem=='ICRS'):
                            print "Reference Sky Coordinate system not J2000 nor ICRS. Quitting."
                            raise

                        (rastring,decstring) = au.rad2radec(rarad,decrad).replace(' ','').split(',')
                        if 'coords' in metadata[band][fname]:
                            metadata[band][fname]['coords'].append([rastring,decstring,refsystem])
                        else:
                            metadata[band][fname]['coords'] = [[rastring,decstring,refsystem]]
                           
                        metadata[band][fname]['coords'].sort()
                        ccss_aux = metadata[band][fname]['coords'][0]
                        ccss_aux_list = [ccss_aux]
                        for ccss_ind in range(1,len(metadata[band][fname]['coords'])):
                            ccss = metadata[band][fname]['coords'][ccss_ind]
                            if not ccss == ccss_aux:
                                ccss_aux_list.append(ccss)
                            ccss_aux = ccss
                            # print "Coordinate mismatch for  field {0}. Filename {1}. Coordinates {2} vs. {3}".\
                            # format(field,filename,ccss_aux,metadata[band][fname]['coords'])
                            # raise
                        metadata[band][fname]['coords'] = ccss_aux_list

            ran_findcont = False
            pipescript = glob.glob("../script/*casa_pipescript.py")
            if len(pipescript) > 0:
                for pscr in pipescript:
                    with open(pscr, 'r') as fh:
                        txt = fh.read()
                    if 'findcont' in txt:
                        ran_findcont = True
                        
            ran_findcont = "ran_findcont" if ran_findcont else "did_not_run_findcont"

            tb.open(filename+"/ANTENNA")
            positions = tb.getcol('POSITION')
            tb.close()
            baseline_lengths = (((positions[None, :, :]-positions.T[:, :, None])**2).sum(axis=1)**0.5)
            maximum_baseline = np.max(baseline_lengths)
            min_bl = np.min(baseline_lengths)
            max_bl = int(np.max(baseline_lengths))

            lb_threshold = {'B3': 750,
                            'B6': 780,
            }
            array_config = ('7M' if max_bl < 100
                            else '12Mshort' if max_bl < lb_threshold[band]
                            else '12Mlong')

            if 'muid_configs' in metadata[band][field]:
                metadata[band][field]['muid_configs'][array_config] = muid
                metadata[band][field]['max_bl'][muid] = max_bl
                metadata[band][field]['min_bl'][muid] = min_bl
            else:
                metadata[band][field]['muid_configs'] = {array_config: muid}
                metadata[band][field]['max_bl'] = {muid: max_bl}
                metadata[band][field]['min_bl'] = {muid: min_bl}

            # Maximum widths to avoid bandwidth smearing
            targetwidth = 0.237 * antenna_diameters[array_config] * bands[band][0] * 1e9 / maximum_baseline
            widths = []
            logprint("Determining smoothing widths for continuum data. targetwidth is {0}".format(targetwidth))
            
            for spw in spws:
                chwid = np.abs(np.mean(msmd.chanwidths(spw)))
                logprint("Channel width for spw {0} is {1}".format(spw,chwid))
                wid = int(targetwidth/chwid)
                if wid <= 0:
                    logprint("Failure at chwid = {0}, wid = {1}.  ".format(chwid, wid))
                    raise ValueError("The channel width is greater than "
                                     "the target line width for spw {0} "
                                     "in ms {1}".format(spw, visfile))
                else:# 
                    # CASA *cannot* handle wid > nchan
                    # This one also insists that there will be at least 3
                    # output channels in all cases. Also, divides nchan.
                    wid = int(msmd.nchan(spw) / 3)
                    divisors = np.array([1])
                    _nn = msmd.nchan(spw)
                    for _pp in [2,3,5]:
                        _prodaux = np.outer(divisors,_pp**np.where(_nn%_pp**np.arange(0,int(ceil(log(_nn,_pp))))==0)[0])
                        divisors = np.sort(np.ndarray.flatten(_prodaux))
                    divisors = divisors[np.where((divisors < wid) & (divisors > _nn/3))]
                    if len(divisors)>0:
                        wid = divisor.max()
                widths.append(wid)

            if 'continuum_widths_in_channels' in metadata[band][field]:
                metadata[band][field]['continuum_widths_in_channels'].append(widths)
            else:
                metadata[band][field]['continuum_widths_in_channels'] = [widths]

            # Custom cont.dat files:
            # <field>.<band>.<array>.cont.dat takes priority; if that exists, it will be used
            # else if
            # <field>.<band>.cont.dat exists, it will be used.
            # we only have 12m and 7m now; everything is otherwise merged
            # (though maybe we'll merge further still)
            arrayname = '12m' if '12M' in array_config else '7m'
            contfile = os.path.join(os.getenv('ALMAIMF_ROOTDIR'),
                                    'contdat',
                                    "{field}.{band}.{array}.cont.dat".format(field=field, band=band,
                                                                             array=arrayname))
            if os.path.exists(contfile):
                logprint("##### Found manually-created cont.dat file {0}".format(contfile))
            else:
                contfile = os.path.join(os.getenv('ALMAIMF_ROOTDIR'),
                                        'contdat',
                                        "{field}.{band}.cont.dat".format(field=field, band=band))
                if os.path.exists(contfile):
                    logprint("##### Found manually-created cont.dat file {0}".format(contfile))
                else:
                    contfile_default = os.path.join(dirpath, '../calibration/cont.dat')
                    logprint("""##### Using default cont.dat file from {0}. 
                        Separating files into fields and copying them into {1}/contdat .""".format(contfile_default,os.getenv('ALMAIMF_ROOTDIR')))
                    #perl stuff to separate cont.dat.
                    os.system(" ".join(['perl', '-lane',
                        '\'if(/Field:\s*(\S+)/){close FF; $aimfdir = $ENV{"ALMAIMF_ROOTDIR"};$field=$1; $nofile=1; if(-e ">$aimfdir/contdat/$1.B6.cont.dat"){$nofile=0}else{open(FF,">$aimfdir/contdat/$1.B6.cont.dat")}}if($field && $nofile){print FF $_}\'', 
                        contfile_default]))
                    contfile = os.path.join(os.getenv('ALMAIMF_ROOTDIR'),
                                    'contdat',
                                    "{field}.{band}.{array}.cont.dat".format(field=field, band=band,
                                                                            array=arrayname))
                    if not os.path.exists(contfile):
                        contfile = os.path.join(os.getenv('ALMAIMF_ROOTDIR'),
                                    'contdat',
                                    "{field}.{band}.cont.dat".format(field=field, band=band,
                                                                            array=arrayname))
 


            if os.path.exists(contfile):
                contdatpath = os.path.realpath(contfile)
                contdat_files[field + band + muid] = contdatpath
                
                cont_channel_selection = parse_contdotdat(contdatpath)
                _, linefracs = contchannels_to_linechannels(cont_channel_selection,
                                                            frqdict,
                                                            return_fractions=True)
                
                if 'cont.dat' in metadata[band][field]:
                    metadata[band][field]['cont.dat'][muid] = contdatpath
                    metadata[band][field]['line_fractions'].append(linefracs)
                else:
                    metadata[band][field]['cont.dat'] = {muid: contdatpath}
                    metadata[band][field]['line_fractions'] = [linefracs]
            else:
                if 'cont.dat' in metadata[band][field]:
                    if muid in metadata[band][field]['cont.dat']:
                        logprint("*** Found DUPLICATE KEY={muid},{max_bl}"
                                 " in cont.dat metadata for band={band} field={field}"
                                 .format(max_bl=max_bl, muid=muid,
                                         band=band, field=field))
                    else:
                        metadata[band][field]['cont.dat'][muid] = 'notfound_' + ran_findcont
                else:
                    metadata[band][field]['cont.dat'] = {muid: 'notfound_' + ran_findcont}
                
                contdat_files[field + band + muid] = 'notfound_' + ran_findcont


            # touch the filename
            with open(os.path.join(dirpath, "{0}_{1}_{2}".format(field, band, array_config)), 'w') as fh:
                fh.write("{0}".format(antnames))
            logprint("Acquired metadata for {0} in {1}_{2}_{3} successfully [Elapsed: {4}]".format(filename, field, band, array_config, time.time() - t0))
                
                
        msmd.close()

with open('metadata.json', 'w') as fh:
    json.dump(metadata, fh)

with open('contdatfiles.json', 'w') as fh:
    json.dump(contdat_files, fh)

logprint("Completed metadata assembly")
logprint("Metadata acquisition took {0} seconds".format(time.time() - t1))
