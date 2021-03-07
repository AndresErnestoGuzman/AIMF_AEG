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
        lines are flmsmdagged out.
    FIELD_ID=<name>
        If this parameter is set, filter out the imaging targets and only split
        fields with this name (e.g., "W43-MM1", "W51-E", etc.).
        Metadata will still be collected for *all* available MSes.
    BAND_TO_IMAGE=B3 or B6
        If this parameter is set, only image the selected band.

The environmental variable ``ALMAIMF_ROOTDIR`` should be set to the directory
containing this file.


Additional Notes
================
USE_SELFCAL_MS is an environmental variable you can set if you want the imaging
to be done using the selfcal.ms file instead of the default continuum MS file.
It is primarily for debug purposes and you shouldn't need it.
"""

onlyDirtyImaging = True
dryRun = False
# PREAMBLE DEFINITIONS
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
            #contimagename = os.path.join(imaging_root, field) +"_"+ band +"_" + arrayname + suffix

            for robust in [0]:
                imname = os.path.join(imaging_root, field) +"_"+ band +"_" + arrayname + suffix+"_robust{0}_dirty".format(robust)
                dirtyImage(visibilities = continuum_visibilities, cell=cell, field=field, band=band, arrayname=arrayname, robust=robust, imaging_root=imaging_root,
                            imsize=imsize, antennae = antennae, phasecenter=phasecenter, suffix=suffix, imaging_parameters = imaging_parameters, dryrun = dryRun)
 
                if onlyDirtyImaging:
                    logprint("Only dirty images this run.")
                    continue
                
                imname = os.path.join(imaging_root, field) +"_"+ band +"_" + arrayname + suffix+"_robust{0}".format(robust)
                image(visibilities=continuum_visibilities, cell=cell, field=field, band=band, arrayname=arrayname, robust=robust, imaging_root=imaging_root,
                        imsize=imsize, antennae = antennae, phasecenter=phasecenter, suffix=suffix, imaging_parameters = imaging_parameters, dryrun = dryRun)

                
