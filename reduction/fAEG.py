
from metadata_tools import determine_imsize, determine_phasecenter, logprint
import re, os, math, copy
from taskinit import msmdtool, iatool, tbtool
from tasks import tclean, exportfits, imstat, imhead, rmtables, split
from utils import validate_mask_path

def create_selfcal_ms_files(mses, arrayname, vis_data, do_bsens=[False]):
    selfcal_mses = []
    tb = tbtool()
    for cms in mses:
        (directory,cmsfilename) = os.path.split(cms)
        scbasename = re.sub('(.+X[^X\.]+).*','\\1',cmsfilename)
        for dbs in do_bsens:
            suffix = "_bsens" if dbs else ""
            scms = "{0}_{1}_selfcal{2}.ms".format(scbasename,arrayname,suffix)
            scms = os.path.join(directory,scms)
            selfcal_mses.append(scms)
            if os.path.exists(scms):
                logprint("Selfcal MS {0} already exists. Using that one.".format(scms),origin='cis_S2_create_selfcal_ms_files')
            else:
                logprint("Selfcal MS {0} does not exists. Splitting a new one.".format(scms),origin='cis_S2_create_selfcal_ms_files')
                tb.open(cms)
                if 'CORRECTED_DATA' in tb.colnames():
                    datacolumn='corrected'
                else:
                    datacolumn='data'
                tb.close()
                split(  vis = cms,
                        outputvis = scms,
                        datacolumn = datacolumn,
                        antenna = vis_data[cms][arrayname+'_antennae'],
                        #spw = spwstr,
                        #width = width,
                        #field = field,
                        )
                logprint("Created new selfcal MS: {0}".format(scms), origin='cis_S2_create_selfcal_ms_files')
    return selfcal_mses


def image(visibilities, cell, field, band, arrayname, robust,imsize, antennae , phasecenter, suffix, imaging_parameters,imaging_root="imaging_results", pbmask=0.25, savemodel='none', datacolumn='corrected',dryrun=False):
    
    logprint("\n*** Starting image function ***",origin='image_function')
    ia = iatool()
    key = "{0}_{1}_{2}_robust{3}".format(field, band,arrayname,robust)
    impars = imaging_parameters[key]
    if not 'maskname' in impars:
        impars['usemask'] = 'pb'
        maskname = ''
        if not 'pblimit' in impars:
            raise KeyError('Either define maskname or set pblimit > 0')

    for iteracion in impars['threshold'].keys():
        impars_thisiter = copy.copy(impars)
        
        for key, val in impars_thisiter.items():
            if isinstance(val, dict):
                impars_thisiter[key] = val[iteracion]

        if impars_thisiter['usemask'] == 'pb'and 'maskname' in impars_thisiter:
            del impars_thisiter['maskname']
        elif impars_thisiter['usemask'] == 'user':
            if 'maskname' in impars_thisiter:
                maskname = validate_mask_path(impars_thisiter['maskname'],
                                              os.getenv('ALMAIMF_ROOTDIR'))
                del impars_thisiter['maskname']
            if 'mask' not in impars_thisiter:
                impars_thisiter['mask'] = maskname

        imname = os.path.join(imaging_root, field) +"_"+ band +"_" + arrayname + suffix+"_robust{0}".format(robust)

        if not os.path.exists(imname+".image.tt0"):
            logprint("Cleaning file {0} for the first time. New image.".format(imname),
                     origin='image_function')
        else:
            logprint("Cleaning already existing file {0}. Will start from this image in tclean.".format(imname),
                     origin='image_function')
        
        
        rmtables(imname+".mask")
        if os.path.exists(imname+".mask"):
            os.system("rm -rf "+imname+".mask")
        remove_repeated_dictionary = [[kk,impars_thisiter.pop(kk,None)] for kk in ['vis','imagename','phasecenter','outframe','veltype','interactive',
                                                    'cell','imsize','antenna','pbcor','savemodel','datacolumn']]
        for keyval in remove_repeated_dictionary:
            if not keyval[1]==None:
                logprint("WARNING: Key {0} with value {1} removed and replaced from 'impars_thisiter' dictionary.".format(keyval[0],keyval[1]),origin='image_function')

        if dryrun:
            logprint("Dry run. tclean("+
            	",".join(["vis={vis}","field={field}","imagename={imagename}","phasecenter={phasecenter}",
            	"outframe={outframe}","veltype={veltype}","interactive={interactive}","cell={cell}",
            	"imsize={imsize}","antenna={antennae}","pbcor={pbcor}","savemodel={savemodel}",
            	"datacolumn={datacolumn}"]).format(vis = visibilities, field = field.encode(), imagename = imname,
                                    phasecenter = phasecenter, outframe = 'LSRK', veltype = 'radio', interactive = False,    
                                    cell = cell, imsize = imsize, antennae = antennae, pbcor = True, savemodel = savemodel,    
                                    datacolumn = datacolumn)+','.join([kk+"="+str(impars_thisiter[kk]) for kk in impars_thisiter.keys()])+")",origin='image_function')
        else:
            tclean(vis = visibilities,
                   field = field.encode(),
                   imagename = imname,
                   phasecenter = phasecenter,
                   outframe = 'LSRK',
                   veltype = 'radio',
                   interactive = False,
                   cell = cell,
                   imsize = imsize,
                   antenna = antennae,
                   pbcor = True,
                   savemodel = savemodel,
                   datacolumn=datacolumn,
                   **impars_thisiter
                  )

            stats = imstat(imagename=imname+".image.tt0",algorithm='hinged-fences',fence=2)
            mad = stats['medabsdevmed'][0]
            rms  = 1.4826* mad
            imhead(imname+".image.tt0",mode='put',hdkey='RMS',hdvalue=rms)

            #os.system("cp -r {0}".image.tt0" {0}.iter{1}".format(imname,iteracion))
            ia.open(imname+".image.tt0")
            ia.sethistory(origin='almaimf_cont_imaging',
               history=["{0}: {1}".format(key, val) for key, val in
                                       impars.items()])
            #ia.sethistory(origin='almaimf_cont_imaging',
            #              history=["git_version: {0}".format(git_version),
            #                       "git_date: {0}".format(git_date)])
            ia.close()

            exportfits(imname+".image.tt0", imname+".image.tt0.fits", overwrite = True)
            exportfits(imname+".image.tt0.pbcor", imname+".image.tt0.pbcor.fits", overwrite = True)
            exportfits(imname+".model.tt0", imname+".model.tt0.fits", overwrite = True)
            exportfits(imname+".residual.tt0", imname+".residual.tt0.fits", overwrite = True)
    logprint("*** Ending image function ***",origin='image_function')

def dirtyImage(visibilities, cell, field, band, arrayname, robust,imsize, antennae , phasecenter, suffix, imaging_parameters,imaging_root="imaging_results", pbmask=0.25,dryrun=False):
    
    logprint("\n*** Starting dirtyImage function ***",origin='dirtyImage_function')
    ia = iatool()
    key = "{0}_{1}_{2}_robust{3}".format(field, band,arrayname,robust)
    impars = imaging_parameters[key]
    dirty_impars = copy.copy(impars)
    if 'maskname' in dirty_impars:
        del dirty_impars['maskname']
    dirty_impars['niter'] = 0
    dirty_impars['usemask'] = 'pb'
    dirty_impars['pbmask'] = pbmask
    if 'threshold' in dirty_impars:
        del dirty_impars['threshold']

    ## Dirty imaging
    imname = os.path.join(imaging_root, field) +"_"+ band +"_" + arrayname + suffix+"_robust{0}_dirty".format(robust)
    if os.path.exists(imname+".image.tt0"):
    	logprint("Skipping dirty image block because {image} already exists. Visibilities={vis}".format(image = imname+
    		".image.tt0",vis = visibilities),origin='dirtyImage_function')
    else:
        logprint("Dirty imaging file {0}".format(imname),
                 origin='almaimf_cont_imaging')

        if dryrun:
            logprint("Dry run. tclean("+
            	",".join(["vis={vis}","field={field}","imagename={imagename}","phasecenter={phasecenter}",
            	"outframe={outframe}","veltype={veltype}","interactive={interactive}","cell={cell}",
            	"imsize={imsize}","antenna={antennae}","pbcor={pbcor}"]).format(vis = visibilities,
            	field = field.encode(), imagename = imname,phasecenter = phasecenter, outframe = 'LSRK',
            	veltype = 'radio', interactive = False, cell = cell, imsize = imsize,
            	antennae = antennae, pbcor = True)+','.join([kk+"="+str(dirty_impars[kk]) for kk in dirty_impars.keys()])+")",origin='dirtyImage_function')
        else:
            tclean(vis = visibilities,
                   field = field.encode(),
                   imagename = imname,
                   #phasecenter = phasecenter,
                   outframe = 'LSRK',
                   veltype = 'radio',
                   #usemask='pb',
                   interactive = False,
                   cell = cell,
                   imsize = imsize,
                   antenna = antennae,
                   pbcor = True,
                   **dirty_impars
                  )

            stats = imstat(imagename=imname+".image.tt0",algorithm='hinged-fences',fence=2)
            mad = stats['medabsdevmed'][0]
            rms  = 1.4826* mad
            imhead(imname+".image.tt0",mode='put',hdkey='RMS',hdvalue=rms)

            ia.open(imname+".residual.tt0")
            ia.sethistory(origin='almaimf_cont_imaging',
                          history=["{0}: {1}".format(key, val) for key, val in
                                   impars.items()])
            #ia.sethistory(origin='almaimf_cont_imaging',
            #              history=["git_version: {0}".format(git_version),
            #                       "git_date: {0}".format(git_date)])
            ia.close()
            exportfits(imname+".image.tt0", imname+".image.tt0.fits", overwrite = True)
    logprint("*** Ending dirtyImage function ***",origin='dirtyImage_function')
        

def imagingOptimalParameters(mses,metadata,exclude_7m = True, only_7m = False):
    
    logprint("\n*** Starting imagingOptimalParameters ***",origin='imagingOptimalParameters_function')
    msmd = msmdtool()
    files_per_field = {}
    vis_image_parameters = {}
    bands = set()
    for ms in mses:
        # get EB name
        ebname =  re.sub('.*(uid[^.]+).*','\\1',os.path.split(ms)[1])
        band = metadata['ebs'][ebname]['band']
        bands = bands | {band}
        vis_image_parameters[ms] = {}
        vis_image_parameters[ms]['ebname'] = ebname

        if not band in files_per_field:
            files_per_field[band] = {}
        fields = metadata['ebs'][ebname]['fields']
        for field in fields:
            if field in files_per_field[band]:
                files_per_field[band][field].append(ms)
            else:
                files_per_field[band][field] = [ms]

        field = fields[0]
        coosys,racen,deccen = determine_phasecenter(ms = ms, field = field)
        # this is center of one of the scecific fields of the visibility. It only approximately 
        # characterizes a single visibility. Possibly OK to caclulate image sizes.
        phasecenter = "{0} {1}deg {2}deg".format(coosys, racen, deccen)

        (dra,ddec,pixscale) = list(determine_imsize(ms=ms, field=field,
                                                    phasecenter=(racen,deccen),
                                                    exclude_7m=exclude_7m,
                                                    only_7m=only_7m,
                                                    spw='all',
                                                    pixfraction_of_fwhm=1/8. if only_7m else 1/4.))
        imsize = [dra, ddec]
        vis_image_parameters[ms]['imsize'] = imsize
        vis_image_parameters[ms]['cellsize'] = ['{0:0.2f}arcsec'.format(pixscale)] * 2
        vis_image_parameters[ms]['cellsize_in_arcsec'] = [math.floor(pixscale*1000 + 0.5)/1000] * 2  # in arcsec
        msmd.open(ms)
        antennae = ",".join([x for x in msmd.antennanames() if 'CM' not in x])
        vis_image_parameters[ms]['12M_antennae'] = antennae
        antennae = ",".join([x for x in msmd.antennanames() if 'CM' in x])
        vis_image_parameters[ms]['7M_antennae'] = antennae
        vis_image_parameters[ms]['7M12M_antennae'] = ""
        msmd.close()


    field_image_parameters = {}
    bands  = files_per_field.keys()
    for band in bands:
        fields = files_per_field[band].keys()
        field_image_parameters[band] = {}
        for field in fields:
            field_image_parameters[band][field] = {}
            field_image_parameters[band][field]['imsize'] = [0,0]
            field_image_parameters[band][field]['cellsize_in_arcsec'] = [1e6,1e6]
            field_image_parameters[band][field]['phasecenter'] = []
            for cvis in files_per_field[band][field]:
                field_image_parameters[band][field]['imsize'] = max(field_image_parameters[band][field]['imsize'],vis_image_parameters[cvis]['imsize'])
                field_image_parameters[band][field]['cellsize_in_arcsec'] = min(field_image_parameters[band][field]['cellsize_in_arcsec'],
                                                                        vis_image_parameters[cvis]['cellsize_in_arcsec'])
                field_image_parameters[band][field]['cellsize'] = ['{0:0.2f}arcsec'.format(field_image_parameters[band][field]['cellsize_in_arcsec'][0])]*2
                coosys,racen,deccen = determine_phasecenter(ms = cvis, field = field)
                phasecenter = "{0} {1}deg {2}deg".format(coosys, racen, deccen)
                field_image_parameters[band][field]['phasecenter'].append(phasecenter)

    logprint("*** Ending imagingOptimalParameters ***",origin='imagingOptimalParameters_function')
    return {'visBased':vis_image_parameters, 'fieldBased':field_image_parameters,'filesPerField': files_per_field}
