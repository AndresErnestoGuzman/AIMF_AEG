from metadata_tools import determine_imsize, determine_phasecenter, logprint, explodeKey
import re, os, math, copy, numpy, shutil
from taskinit import msmdtool, casalog, qatool, tbtool, mstool, iatool
from tasks import tclean, exportfits, imstat, imhead, rmtables, split, delmod, imregrid
from utils import validate_mask_path

def create_clean_model(cubeimagename, contimagename, imaging_results_path='.', contmodel_path=None,nott1=False):
    
    ia = iatool()
    contimagename = re.sub(r"\.image\.tt[0-1](\.fits)?$","",contimagename)
    cubeimagename = re.sub(r"\.image(\.fits)?$","",cubeimagename)

    #results_path = "./imaging_results/"  # imaging_results frmo the pipeline
    #contmodel_path = "./imaging_results_test_casatools/"  #Path with input and temporary continuum models
    if contmodel_path is None:
        contmodel_path = imaging_results_path

    # Create continuum_model.image.tt0 and .tt1 regridded to the cube spatial frame, but still with 1 spectral pix
    
    tt0name = "{contmodelpath}/{contimagename}.model.tt0".format(contimagename=contimagename, contmodelpath=contmodel_path)
    tt1name = "{contmodelpath}/{contimagename}.model.tt1".format(contimagename=contimagename, contmodelpath=contmodel_path)
    if not os.path.exists(tt0name):
        raise IOError("Continuum startmodel file {0} does not exist.".format(tt0name))
    if not os.path.exists(tt1name):
        logprint("Only tt0 for continuum model.",origin="create_clean_model")
        nott1 = True

    temp_dict_cont_tt0 = imregrid(imagename=tt0name, template="get")
    # not needed temp_dict_cont_tt1 = imregrid(imagename=tt1name, template="get")

    # image has to exist, model should not exist!
    # (if you ran tclean with niter=0, no model is created)
    cubeinimagepath = ("{results_path}/{cubeimagename}.image"
                       .format(results_path=imaging_results_path,
                               cubeimagename=cubeimagename))
    cubeoutmodelpath = ("{results_path}/{cubeimagename}.contcube.model"
                        .format(results_path=imaging_results_path,
                                cubeimagename=cubeimagename))

    if not os.path.exists(cubeinimagepath):
        raise IOError("Cube file {0} does not exist.".format(cubeinimagepath))

    temp_dict_line = imregrid(imagename=cubeinimagepath, template="get")
    temp_dict_line['shap'][-1] = 1
    temp_dict_line['csys']['spectral2'] = temp_dict_cont_tt0['csys']['spectral2']
    temp_dict_line['csys']['worldreplace2'] = temp_dict_cont_tt0['csys']['worldreplace2']

    tt0model = ("{contmodel_path}/{cubeimagename}_continuum_model.image.tt0"
                .format(contmodel_path=contmodel_path,
                        cubeimagename=cubeimagename))
    tt1model = ("{contmodel_path}/{cubeimagename}_continuum_model.image.tt1"
                .format(contmodel_path=contmodel_path,
                        cubeimagename=cubeimagename))

    imregrid(imagename=tt0name, output=tt0model, template=temp_dict_line, overwrite=True)
    imregrid(imagename=tt1name, output=tt1model, template=temp_dict_line, overwrite=True)

    # Use CASA tools to create a model cube from the continuum model
    if os.path.exists(cubeoutmodelpath):
        shutil.rmtree(cubeoutmodelpath)
    shutil.copytree(cubeinimagepath, cubeoutmodelpath)

    dict_line = imregrid(imagename=cubeoutmodelpath, template="get")
    line_im = ia.newimagefromfile(cubeoutmodelpath)
    #print(line_im.shape())

    tt0_im = ia.newimagefromfile(tt0model)
    #print(tt0_im.shape())
    tt0_pixvalues = tt0_im.getchunk()
    tt0_im.close()

    tt1_im = ia.newimagefromfile(tt1model)
    #print(tt1_im.shape())
    tt1_pixvalues = tt1_im.getchunk()
    tt1_im.close()

    # From Eq. 2 of Rau & Cornwell (2011)
    # temp_dict_cont_tt0['csys']['spectral2']['wcs']
    # dict_line['csys']['spectral2']['wcs']
    # dnu_plane: dnu with respect to cube reference freq.
    # dnu: dnu with respect to tt0 continuum reference
    for plane in range(dict_line['shap'][-1]):
        if (not plane%50) or plane == dict_line['shap'][-1]:
            logprint('Calculating continuum model for plane {}'.format(plane))
        dnu_plane = (plane - dict_line['csys']['spectral2']['wcs']['crpix'])*dict_line['csys']['spectral2']['wcs']['cdelt']
        nu_plane = dict_line['csys']['spectral2']['wcs']['crval'] + dnu_plane
        #print(dnu_plane, nu_plane)
        if nott1:
            plane_pixvalues = tt0_pixvalues
        else:
            factor = (nu_plane - temp_dict_cont_tt0['csys']['spectral2']['wcs']['crval'])/temp_dict_cont_tt0['csys']['spectral2']['wcs']['crval']
            #print(factor)
            plane_pixvalues = tt0_pixvalues + factor*tt1_pixvalues
        
        blc = [0, 0, 0, plane]
        #trc = [line_im.shape()[0]-1, line_im.shape()[1]-1, 0, plane]
        line_im.putchunk(plane_pixvalues, blc=blc, replicate=False)
    
    line_im.close()

    return cubeoutmodelpath


def line_imaging_dict(chanwidth=None,fullwidth=None,restfreq=None,start=None,nchan=-1,cubeasymmetry=0.5,vlsrkms=0):
    qa = qatool()
    dd = {'nchan': nchan}
    if chanwidth:
        dd['width'] = chanwidth
        chanWidth = qa.quantity(chanwidth)
    if restfreq:
        dd['restfreq'] = restfreq
        restFreq = qa.quantity(restfreq)
    if start:
        dd['start'] = start
        start = qa.quantity(start)
    if fullwidth:
        fullWidth = qa.quantity(fullwidth)
        if not qa.getunit(fullWidth) == 'km/s':
            raise Exception("Units!")
        if restfreq:
            start = qa.quantity((-1.0+cubeasymmetry)*qa.getvalue(fullWidth)+vlsrkms,'km/s')
            dd['start'] = qa.tos(start)            
            if chanwidth:
                fullWidth = qa.convert(fullWidth,qa.getunit(chanwidth))
                dd['nchan'] = int(qa.getvalue(qa.div(fullWidth,chanWidth)))
            return dd

    return dd


def rms_from_mad(imagename):
        stats = imstat(imagename=imagename,algorithm='hinged-fences',fence=2)
        mad = stats['medabsdevmed'][0]
        rms  = 1.4826* mad
        return(rms)


def selfcal_name(ms, arrayname, field, band, bsens=False):
    (directory,msfilename) = os.path.split(ms)
    suffix = "_bsens" if bsens else ""
    scbasename = re.sub('(.+X[^X\.]+).*','\\1',msfilename)
    scms = "{uid}_{fi}_{band}_{arr}_selfcal{suffix}.ms".format(uid=scbasename,arr=arrayname,suffix=suffix,fi=field,band=band)
    scms = os.path.join(directory,scms)
    return scms

class UniversalSet(set):
    def __and__(self, other):
        return other

    def __rand__(self, other):
        return other

def sc_image(visibilities, cell, field, band, arrayname, robust,imsize, antennae , phasecenter, suffix, imaging_parameters ,imagename=None,imaging_root="imaging_results", pbmask=0.25, savemodel='none', datacolumn='corrected',dryrun=False):
    
    logprint("\n*** Starting image function ***",origin='sc_image_function')
    ia = iatool()
    key = "{0}_{1}_{2}_robust{3}".format(field, band,arrayname,robust)
    if key in imaging_parameters:
        impars = imaging_parameters[key]
    else:
        impars = imaging_parameters

    if not 'maskname' in impars:
        impars['usemask'] = 'pb'
        maskname = ''
        if not 'pblimit' in impars:
            raise Exception('Either define maskname or set pblimit > 0')

    if impars['usemask'] == 'pb'and 'maskname' in impars:
        del impars['maskname']
    elif impars['usemask'] == 'user':
        if 'maskname' in impars:
            maskname = validate_mask_path(impars['maskname'],
                                          os.getenv('ALMAIMF_ROOTDIR'))
            del impars['maskname']
        if 'mask' not in impars:
            impars['mask'] = maskname

    if imagename==None:
        imname = os.path.join(imaging_root, field +"_"+ band +"_" + arrayname + "_" + suffix+"robust{0}".format(robust))
    else:
        imname = imagename

    if not os.path.exists(imname+".image.tt0"):
        logprint("Cleaning file {0} for the first time. New image.".format(imname),
                 origin='image_function')
    else:
        logprint("Cleaning already existing file {0}. Will start from this image in tclean.".format(imname),
                 origin='image_function')
    
    
    rmtables(imname+".mask")
    if os.path.exists(imname+".mask"):
        os.system("rm -rf "+imname+".mask")
    remove_repeated_dictionary = [[kk,impars.pop(kk,None)] for kk in ['vis','imagename','phasecenter','outframe','veltype','interactive',
                                                'cell','imsize','antenna','pbcor','savemodel','datacolumn']]
    for keyval in remove_repeated_dictionary:
        if not keyval[1]==None:
            logprint("WARNING: Key {0} with value {1} removed and replaced from 'impars' dictionary.".format(keyval[0],keyval[1]),origin='sc_image_function')

    if dryrun:
        logprint("Dry run. tclean("+
            ",".join(["vis={vis}","field={field}","imagename={imagename}","phasecenter={phasecenter}",
            "outframe={outframe}","veltype={veltype}","interactive={interactive}","cell={cell}",
            "imsize={imsize}","antenna={antennae}","pbcor={pbcor}","savemodel={savemodel}",
            "datacolumn={datacolumn}"]).format(vis = visibilities, field = field.encode(), imagename = imname,
                                phasecenter = phasecenter, outframe = 'LSRK', veltype = 'radio', interactive = False,    
                                cell = cell, imsize = imsize, antennae = antennae, pbcor = True, savemodel = savemodel,    
                                datacolumn = datacolumn)+','.join([kk+"="+str(impars[kk]) for kk in impars.keys()])+")",origin='sc_image_function')
        for visi in visibilities:
            if not os.path.exists(visi):
                logprint("{0} does not exists!. Terminating dry run".format(visi),origin='sc_image_function')
                raise Exception("{0} does not exists!. Terminating dry run".format(visi))        
    else:
        tcleanargs = getattr(tclean,'parameters').keys()
        difargs = set(impars.keys()) - set(tcleanargs)
        for difarg in difargs:
            del impars[difarg]
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
               **impars
              )

        stats = imstat(imagename=imname+".image.tt0",algorithm='hinged-fences',fence=2)
        mad = stats['medabsdevmed'][0]
        rms  = 1.4826* mad
        imhead(imname+".image.tt0",mode='put',hdkey='RMS',hdvalue=rms)

        ia.open(imname+".image.tt0")
        ia.sethistory(origin='image_function',
           history=["{0}: {1}".format(key, val) for key, val in impars.items()])
        ia.close()

        exportfits(imname+".image.tt0", imname+".image.tt0.fits", overwrite = True)
        exportfits(imname+".image.tt0.pbcor", imname+".image.tt0.pbcor.fits", overwrite = True)
        exportfits(imname+".model.tt0", imname+".model.tt0.fits", overwrite = True)
        exportfits(imname+".residual.tt0", imname+".residual.tt0.fits", overwrite = True)
    logprint("*** Ending image function ***",origin='sc_image_function')

def sc_iter_image(iteracion,arrayname,selfcal_files_per_field,do_bsens, field_image_parameters,vis_image_parameters,continuum_files_per_field,selfcal_imaging_pars,imaging_root,dryRun ,savemodel,datacolumn,fields=UniversalSet(), bands=UniversalSet()):
    
    imagenames = []
    bands = set(selfcal_files_per_field.keys()) & bands
    for band in bands:
        fields = set(selfcal_files_per_field[band].keys()) & fields
        if not os.getenv('FIELD_ID') == None:
            fields = fields & set(os.getenv('FIELD_ID').split())
        logprint("Doing band {band} images for fields {fields}".format(band=band,fields=", ".join(list(fields))), origin='sc_iter_image')
        for field in fields:
            dbs = do_bsens
            suffix = "_bsens" if dbs else ""
            phasecenter = field_image_parameters[band][field]['phasecenter'][0]
            antennae = [vis_image_parameters[vis][arrayname+'_antennae'] for vis in continuum_files_per_field[band][field]]

            selfcal_visibilities = [x.replace("selfcal.ms","selfcal"+suffix+".ms") for x in selfcal_files_per_field[band][field]]
            
            imsize = field_image_parameters[band][field]['imsize']
            cell = field_image_parameters[band][field]['cellsize']

            for robust in [0]:
                key = "{0}_{1}_{2}_robust{3}".format(field, band,arrayname,robust)
                keybs = "{0}_{1}_{2}_bsens_robust{3}".format(field, band,arrayname,robust)
                if dbs and (keybs in selfcal_imaging_pars):
                    key = keybs

                selfcal_imaging_pars_thisiter = select_iter(selfcal_imaging_pars[key],iteracion)

                imname = os.path.join(imaging_root, field) +"_"+ band +"_" + arrayname + suffix+"_robust{0}_selfcal{1}".format(robust,iteracion)
                imagenames.append(imname)
                if iteracion==0 or os.path.exists(imname+".image.tt0"):
                    selfcal_imaging_pars_thisiter['startmodel'] = ''
                    for vv in selfcal_visibilities:
                        delmod(vis=vv)
                else:
                    selfcal_imaging_pars_thisiter['startmodel'] = [imname.replace("selfcal"+str(iteracion),"selfcal"+str(iteracion-1))+".model.tt"+str(ttn) for ttn in [0,1]]

                if os.path.exists(imname+".image.tt0"):
                    logprint("Image {0} already exists. Starting from this one.".format(os.path.split(imname)[1]),origin='sc_iter_image')
                
                logprint("Doing selfcal. {0} image.".format(os.path.split(imname)[1]),origin='sc_iter_image')
                sc_image(visibilities=selfcal_visibilities, cell=cell, field=field, band=band, arrayname=arrayname, robust=robust, imaging_root=imaging_root,
                        imsize=imsize, antennae = antennae, phasecenter=phasecenter, suffix=suffix, imaging_parameters = selfcal_imaging_pars_thisiter, dryrun = dryRun,
                        savemodel=savemodel,datacolumn=datacolumn,imagename=imname)

                if (not dryRun) and iteracion==0:
                    if all(check_model_column(vis=selfcal_visibilities)):
                        logprint("Model column populated from previous selfcal step.",origin='sc_iter_image_check_model_column')
    return(imagenames)

def select_iter(dictionary,i):
    rr = {}
    for k in dictionary.keys():
        if isinstance(dictionary[k],dict):
            if i in dictionary[k]:
                rr[k] = dictionary[k][i]
            elif str(i) in dictionary[k]:
                rr[k] = dictionary[k][str(i)]
        else:
            rr[k] = dictionary[k]
    return rr

def check_model_column(vis):
    
    if isinstance(vis,str):
        vis = [vis]
    elif not isinstance(vis,list):
        logprint("Visibility type {0} not recognized.".format(vis),origin=check_model_column)
        raise Exception

    rr =[]
    for visi in vis:
        ms = mstool()
        ms.open(visi)
        model_data = ms.getdata(['MODEL_PHASE'])
        ms.close()
        if 'model_phase' not in model_data:
            logprint("Error encountered: model column of {0} was not populated!".format(os.path.split(visi)[1]),
                     origin='check_model_column')
            raise Exception
        elif numpy.all(model_data['model_phase'] == 0):
            logprint("Model column of {0} is zero!".format(os.path.split(visi)[1]),
                     origin='check_model_column')
            raise Exception
        else:
            logprint("Model column of {0} is populated.".format(os.path.split(visi)[1]),
                     origin='check_model_column')
        rr.append(True)
    return rr

def create_selfcal_ms_files(mses, arrayname, vis_data,eb_data, do_bsens=[False]):
    selfcal_mses = []
    tb = tbtool()
    for cms in mses:
        (directory,msfilename) = os.path.split(cms)
        ebname = re.sub('.*(uid.+X[^X\.]+).*','\\1',msfilename)
        fields = eb_data[ebname]['fields']
        band = eb_data[ebname]['band']
        for dbs in do_bsens:
            suffix  = '_bsens' if dbs else ''
            visi = cms.replace("split.cal.cont","split.cal"+suffix+".cont")
            for field in fields:
                scms = selfcal_name(ms = visi, arrayname = arrayname, field = field, band = band, bsens = dbs)
                
                # for key in selfcal_mses_data[scms].keys():
                #     if 'antennae' in key and not key==arrayname+'_antennae':
                #         del  selfcal_mses_data[scms][key]

                if os.path.exists(scms):
                    logprint("Selfcal MS {0} already exists. Using that one.".format(scms),origin='create_selfcal_ms_files')
                else:
                    logprint("Selfcal MS {0} does not exists. Splitting a new one.".format(scms),origin='create_selfcal_ms_files')
                    tb.open(visi)
                    if 'CORRECTED_DATA' in tb.colnames():
                        datacolumn='corrected'
                    else:
                        datacolumn='data'
                    tb.close()
                    split(  vis = visi,
                            outputvis = scms,
                            datacolumn = datacolumn,
                            antenna = vis_data[cms][arrayname+'_antennae'],
                            #spw = spwstr,
                            #width = width,
                            field = field,
                            )
                    logprint("Created new selfcal MS: {0} from {1}".format(scms, visi), origin='create_selfcal_ms_files')
                selfcal_mses.append(scms)
    return selfcal_mses

def image(visibilities, cell, field, band, arrayname, robust,imsize, antennae , phasecenter, suffix, imaging_parameters ,imagename=None,imaging_root="imaging_results", pbmask=0.25, savemodel='none', datacolumn='corrected',dryrun=False):
    
    logprint("\n*** Starting image function ***",origin='image_function')
    ia = iatool()
    key = "{0}_{1}_{2}_robust{3}".format(field, band,arrayname,robust)
    if key in imaging_parameters:
        impars = imaging_parameters[key]
    else:
        impars = imaging_parameters

    if not 'maskname' in impars:
        impars['usemask'] = 'pb'
        maskname = ''
        if not 'pblimit' in impars:
            raise Exception('Either define maskname or set pblimit > 0')

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

        if imagename==None:
            imname = os.path.join(imaging_root, field) +"_"+ band +"_" + arrayname + suffix+"_robust{0}".format(robust)
        else:
            imname = imagename

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
            tcleanargs = getattr(tclean,'parameters').keys()
            difargs = set(impars_thisiter.keys()) - set(tcleanargs)
            for difarg in difargs:
                del impars_thisiter[difarg]
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
            ia.sethistory(origin='image_function',
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

def dirtyImage(visibilities, cell, field, band, arrayname, robust,imsize, antennae , phasecenter, suffix, imaging_parameters,imagename=None,imaging_root="imaging_results", pbmask=0.25,dryrun=False):
    
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
    if imagename==None:
        imname = os.path.join(imaging_root, field) +"_"+ band +"_" + arrayname + suffix+"_robust{0}_dirty".format(robust)
    else:
        imname = imagename

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
            tcleanargs = getattr(tclean,'parameters').keys()
            difargs = set(dirty_impars.keys()) - set(tcleanargs)
            for difarg in difargs:
                del dirty_impars[difarg]
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
            ia.sethistory(origin='dirtyImage_function',
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


