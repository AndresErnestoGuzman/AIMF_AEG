"""
Imaging parameters for all continuum imaging work.

The first variable is used by ``continuum_imaging.py``.  It specifies the
parameters to be used for the first-pass imaging before any self-calibration is
done.  Please add your source name to the 'field' section.
DO NOT modify the default imaging parameters; if there are any that are
unsuitable for your target field, add them to the
``imaging_parameters_nondefault`` keyword, following the naming scheme laid out
there.

If you want to specify a different set of imaging parameters, you can do so by
passing a dictionary instead of a single number.  For example, instead of
    threshold: '1mJy'
you can use
    threshold: {0: '2mJy', 1:'1mJy', 2:'0.75mJy'}
The requirements are:
    (1) you must have a zero entry (which is used by continuum_imaging.py)
    (2) you must have the same number of entries in each dictionary as there
    are in the calibration parameters list below

The self-calibration parameters are specified in ``selfcal_pars``.  The default is to
do 4 iterations of phase-only self calibration.  If you would like to add additional
steps, you can add them by adding new entries to the self-calibration parameter
dictionary for your source following the template laid out below.


You can copy any set of parameters and add `_bsens` to the end of the name to
have it use special parameters only for the bsens imaging.  For example:
    'G343.1261-00.0623_B3_12M_robust0': {...},
    'G343.1261-00.0623_B3_12M_robust0_bsens: {...},'
would be the parameters used for non-bsens and for bsens data, respectively.
If you have ONLY a non-bsens parameter key (you do not have a _bsens set
of parameters), the bsens selfcal & imaging will use the same as the non-bsens.

CONTRIBUTOR NOTE:
    This file is to be formatted with python's "black" formatter:

        black -t py27 -l 120 imaging_parameters.py
"""

import os, sys, glob, subprocess, re, copy, json
from metadata_tools import json_load_byteified, explodeKey

if os.path.exists('metadata.json'):
    with open('metadata.json', 'r') as fh:
        metadata = json_load_byteified(fh)
else:
    raise FileError("metadata.json not found!")


allfields = set()
for band in set(metadata.keys()).difference({'fields','ebs'}):
    allfields = allfields | set(str(x) for x in metadata[band])
# if 'FIELD_ID' in os.environ:
#     allfields = allfields & set(os.getenv('FIELD_ID').split())
print allfields

# set up global defaults. Mainly for 
imaging_parameters = {
    "{0}_{1}_{2}_robust{3}".format(field, band, array, robust): {
        #"threshold": "1mJy",  # RMS ~0.5-0.6 mJy
        "pblimit": 0.1,
        "niter": 1000,
        "robust": robust,
        "weighting": "briggs",
        "scales": [0, 3, 9],
        #"gridder": "mosaic",
        "specmode": "mfs",
        "deconvolver": "mtmfs",
        "usemask": "pb",
        "nterms": 2,
    }
    for field in allfields
    for band in ("B3", "B6")
    for array in ("12M", "7M12M", "7M")
    for robust in [0]
}


# added for 7M only data: higher threshold
for key in imaging_parameters:
    dd = explodeKey(key)
    if dd['array'] == "7M":
        imaging_parameters[key]["threshold"] = {0:"5mJy"}
        imaging_parameters[key]["scales"] = [0, 3, 9, 27]

imaging_parameters_nondefault = {}

if True:
    mascaras = glob.glob("AIMF_AEG/reduction/clean_regions/*M_robust0.image.tt0_4.0sigma_*Jy")
    #re.sub('.*(uid[^.]+).*','\\1',os.path.split(continuum_ms)[1])
    for mm in mascaras:
        (path,filename) = os.path.split(mm)
        match = re.match(r"([^_]+)_(B[1-9]+)_([127]+M)_.*(robust[^\._]+).*\.image\.tt0",filename)
        basename = match.group(0)
        (field, band, array, robust) = match.group(1,2,3,4)
        key = '_'.join(match.group(1,2,3,4))
        imaging_parameters_nondefault[key] = {}

        nsigma = 4
        maskname = glob.glob("AIMF_AEG/reduction/clean_regions/"+basename+"_"+('{0:.1f}'.format(nsigma))+"sigma*Jy")
        if len(maskname)>0:
            maskname = maskname[0].replace('.fits','')
            (path,maskname) = os.path.split(maskname)
            auxiliar = maskname.split('_')
            thresholdmJy = auxiliar.pop()
            thresholdmJy = float(thresholdmJy.replace("mJy",""))
            sigma = thresholdmJy*1e-3/nsigma
            imaging_parameters_nondefault[key]['sigma'] = sigma
            imaging_parameters_nondefault[key]['threshold'] = {0: "{0:.3f}mJy".format(2.0e3*sigma), 1:"{0:.3f}mJy".format(2.0e3*sigma)}
            imaging_parameters_nondefault[key]['niter']     = {0: 6000,  1:12000}
            imaging_parameters_nondefault[key]['usemask']   = {0:'user', 1:'pb'}
            imaging_parameters_nondefault[key]['maskname']  = {0: maskname,1:''}
            imaging_parameters_nondefault[key]['pbmask']  = {0: 0.0,1:0.25}
        elif os.path.exists('imaging_results/'+basename+'.fits'):
            sigma = float(subprocess.check_output(['gethead','imaging_results/'+basename+'.fits','RMS']))
            #sigma = 1
            imaging_parameters_nondefault[key]['sigma'] = sigma
            imaging_parameters_nondefault[key]['threshold'] = {0: "{0:.3f}mJy".format(2.0e3*sigma)}
            imaging_parameters_nondefault[key]['niter']     = {0: 6000}
            imaging_parameters_nondefault[key]['usemask']   = {0:'pb'}
            imaging_parameters_nondefault[key]['maskname']  = {0: ''}
            imaging_parameters_nondefault[key]['pbmask']  = {0: 0.0,1:0.25}
        else:
            print "Cannot determine sigma nor threshold for {0} sir.\nFinishing.".format(basename)
            raise

for key in imaging_parameters_nondefault:
    if "bsens" in key:
        check_key = "_".join(key.split("_")[:-1])
        assert check_key in imaging_parameters, "key {0} not in imaging_parameters!".format(check_key)
        imaging_parameters[key] = copy.deepcopy(imaging_parameters[check_key])
    else:
        assert key in imaging_parameters, "key {0} was not in imaging_parameters".format(key)
    imaging_parameters[key].update(imaging_parameters_nondefault[key])

# Default sigma calculation for threshold
if False:
    images_for_sigma_estimation = glob.glob("imaging_results/*M_robust0.image.tt0.fits")
    for i in images_for_sigma_estimation:
        (path,filename) = os.path.split(i)
        auxiliar = filename.split('_')
        field = auxiliar.pop(0)
        band = auxiliar.pop(0)
        array = auxiliar.pop(0)
        robust_value = re.sub('robust([^\.]+).*','\\1',auxiliar.pop(0))
        key = "{0}_{1}_{2}_robust{3}".format(field, band, array, robust_value)
        if not 'threshold' in imaging_parameters[key]:
            rms = float(subprocess.check_output(['gethead','RMS',i]))
            # 2 sigma threshold
            imaging_parameters[key]['threshold'] = {0:"{0:.2f}mJy".format(rms*1000*2)}


sc_ip = copy.deepcopy(imaging_parameters)
images = glob.glob("AIMF_AEG/reduction/clean_regions/*M_robust0.image.tt0_4.0sigma_*Jy")
for mm in images:
    (path,filename) = os.path.split(mm)
    match = re.match(r"([^_]+)_(B[1-9]+)_([127]+M)_.*(robust[^\._]+).*\.image\.tt0",filename)
    basename = match.group(0)
    (field, band, array, robust) = match.group(1,2,3,4)
    key = '_'.join(match.group(1,2,3,4))
    imaging_parameters_nondefault[key] = {}

    masknames = glob.glob("AIMF_AEG/reduction/clean_regions/"+basename+"*sigma*Jy")
    masknames.sort(reverse=True)
    sc_ip[key]['maskname'] = {}
    sc_ip[key]['niter'] = {}
    sc_ip[key]['threshold'] = {}
    sc_ip[key]['usemask'] = {}
    sc_ip[key]['pbmask'] = {}
    i = 0
    for mm in masknames:
        match = re.match(r".*_([0-9\.]+)sigma.*",mm)
        nsigma = match.group(1)
        sc_ip[key]['maskname'][i] = mm
        sc_ip[key]['niter'][i] = 2**(i+5)
        # clean up to half the mask limit
        sc_ip[key]['threshold'][i] = "{0:.3f}mJy".format(float(nsigma)*sc_ip[key]['sigma']*1000./2) 
        sc_ip[key]['usemask'][i] = 'user'
        sc_ip[key]['pbmask'][i] = 0.0
        i += 1
    
    sc_ip[key]['maskname'][i] = ''
    sc_ip[key]['niter'][i] = 7000
    # clean up to half the mask limit
    sc_ip[key]['threshold'][i] = "{0:.3f}mJy".format(sc_ip[key]['sigma']*2.e3)
    sc_ip[key]['usemask'][i] = 'pb'
    sc_ip[key]['pbmask'][i] = 0.25

    i += 1
    sc_ip[key]['maskname'][i] = ''
    sc_ip[key]['niter'][i] = 7000
    # clean up to half the mask limit
    sc_ip[key]['threshold'][i] = "{0:.3f}mJy".format(sc_ip[key]['sigma']*2.e3)
    sc_ip[key]['usemask'][i] = 'pb'
    sc_ip[key]['pbmask'][i] = 0.25


selfcal_imaging_pars = sc_ip
"""
Self-calibration parameters are defined here
"""

default_selfcal_pars = {
        "solint": "inf",
        "gaintype": "T",
        "solnorm": True,
        # 'combine': 'spw', # consider combining across spw bounds
        "calmode": "p",}

sc_p = {}
for key in sc_ip.keys():
    if '7M' in key or not 'B6' in key:
        continue
    sc_p[key] = {}
    for i in sorted(sc_ip[key]['threshold'].keys())[1:]: # Starting in 1.
        sc_p[key][i] = copy.deepcopy(default_selfcal_pars)
        if i==1:
            sc_p[key][i] = {"calmode": "p", "gaintype": "G", "solint": "inf", "solnorm": True}
        if i>1:
            sc_p[key][i]['solint'] =  "{0}s".format(10+10*(max(sc_ip[key]['threshold'].keys())-i))
        if i == max(sc_ip[key]['threshold'].keys()):
            sc_p[key][i] = {"calmode": "ap", "gaintype": "T", "solint": "inf", "solnorm": True}

selfcal_pars_custom = {
#     "G343.1261-00.0623_B6_12M_robust-2": {
#         1: {"calmode": "p", "gaintype": "T", "solint": "inf", "solnorm": True},
#         2: {"calmode": "p", "gaintype": "T", "solint": "inf", "solnorm": True},
#         3: {"calmode": "p", "gaintype": "T", "solint": "inf", "solnorm": True},
#         4: {"calmode": "p", "gaintype": "T", "solint": "inf", "solnorm": True},
#     },
 
}


selfcal_pars = sc_p.copy()
for key in selfcal_pars_custom:
    for iternum in selfcal_pars_custom[key]:
        if iternum in selfcal_pars[key]:
            selfcal_pars[key][iternum].update(selfcal_pars_custom[key][iternum])
        else:
            selfcal_pars[key][iternum] = selfcal_pars_custom[key][iternum]

## LINE STUFF

# line_imaging_parameters_default = {
#     "{0}_{1}_{2}_robust{3}{4}".format(field, band, array, robust, contsub): {
#         "niter": 5000000,
#         "threshold": "5sigma",
#         "robust": robust,
#         "weighting": "briggs",
#         "deconvolver": "hogbom",
#         # "scales": [0, 3, 9, 27, 81],
#         # "nterms": 1,
#         "gridder": "mosaic",
#         "specmode": "cube",
#         "outframe": "LSRK",
#         "veltype": "radio",
#         "usemask": "pb",
#         "pblimit": 0.05,
#         "pbmask": 0.1,
#         "perchanweightdensity": False,
#         "interactive": False,
#         "mask_out_endchannels": 2,
#     }
#     for field in allfields
#     for band in ("B3", "B6")
#     for array in ("12M", "7M12M", "7M")
#     # for robust in (0,)
#     for robust in [0]
#     for contsub in ("", "_contsub")
# }

# line_imaging_parameters = copy.deepcopy(line_imaging_parameters_default)

# line_imaging_parameters_custom = {
 
#     "G343.1261-00.0623_B3_12M_robust0": {
#         "threshold": "6mJy",
#         # "startmodel": "G343.1261-00.0623_B3_uid___A001_X1296_X1e9_continuum_merged_12M_robust0_selfcal5_finaliter",
#         # UF machine has 1e9 instead of 1e5 as of 10/14/2020 - this may change
#         "startmodel": "G343.1261-00.0623_B3_uid___A001_X1296_X1e5_continuum_merged_12M_robust0_selfcal7_finaliter",
#     },
#     "G343.1261-00.0623_B3_12M_robust0_contsub": {  # More restrictive params for bright emission
#         "usemask": "auto-multithresh",
#         "noisethreshold": 4.2,
#         "sidelobethreshold": 2.5,
#         "negativethreshold": 15.0,
#         "minbeamfrac": 0.3,
#         "dogrowprune": True,
#         "threshold": "6mJy",
#     },
#     "G343.1261-00.0623_B6_12M_robust0": {
#         "threshold": "6mJy",
#         # "startmodel": "G343.1261-00.0623_B6_uid___A001_X1296_X1df_continuum_merged_12M_robust0_selfcal5_finaliter",
#         # UF machine has 1db instead of 1df as of 10/16/2020 - this may change
#         "startmodel": "G343.1261-00.0623_B6_uid___A001_X1296_X1db_continuum_merged_12M_robust0_selfcal5_finaliter",
#     }
# }


# for key in line_imaging_parameters_custom:
#     if key in line_imaging_parameters:
#         line_imaging_parameters[key].update(line_imaging_parameters_custom[key])
#     else:
#         line_imaging_parameters[key] = line_imaging_parameters_custom[key]

# default_lines = {
#     "n2hp": "93.173700GHz",
#     "sio": "217.104984GHz",
#     "h2co303": "218.222195GHz",
#     "12co": "230.538GHz",
#     "13cs": "231.22068520GHz",
#     "h30a": "231.900928GHz",
#     "h41a": "92.034434GHz",
#     "c18o": "219.560358GHz",
#     "ch3cn": "92.26144GHz",
#     "ch3cch": "102.547983GHz",
# }
# field_vlsr = {
#     "G343.1261-00.0623": "-2km/s",
# }
# # line parameters are converted by line_imaging.py into tclean parameters
# line_parameters_default = {
#     field: {
#         line: {"restfreq": freq, "vlsr": field_vlsr[field], "cubewidth": "50km/s"}
#         for line, freq in default_lines.items()
#     }
#     for field in allfields
# }
# for field in allfields:
#     line_parameters_default[field]["12co"]["cubewidth"] = "150km/s"
#     line_parameters_default[field]["ch3cn"]["cubewidth"] = "150km/s"  # is 150 wide enough?
# line_parameters = copy.deepcopy(line_parameters_default)

# line_parameters_custom = {
 
#     "G343.1261-00.0623": {
#         "12co": {"cubewidth": "150km/s"},
#         "ch3cn": {"cubewidth": "150km/s"},
#         "h41a": {"cubewidth": "120km/s", "vlsr": "0km/s", "width": "2km/s"},
#         "n2hp": {"cubewidth": "60km/s"},
#     }
# }

# for field in line_parameters_custom:
#     for line in line_parameters_custom[field]:
#         line_parameters[field][line].update(line_parameters_custom[field][line])
