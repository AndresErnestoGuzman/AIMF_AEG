from metadata_tools import determine_imsize

def imagingOptimalParameters(mses,metadata):

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
        phasecenter = "{0} {1}deg {2}deg".format(coosys, racen, deccen)
        vis_image_parameters[ms]['phasecenter'] = phasecenter
        (dra,ddec,pixscale) = list(determine_imsize(ms=ms, field=field,
                                                    phasecenter=(racen,deccen),
                                                    exclude_7m=exclude_7m,
                                                    only_7m=only_7m,
                                                    spw='all',
                                                    pixfraction_of_fwhm=1/8. if only_7m else 1/4.))
        imsize = [dra, ddec]
        vis_image_parameters[ms]['imsize'] = imsize
        vis_image_parameters[ms]['cellsize'] = ['{0:0.2f}arcsec'.format(pixscale)] * 2
        vis_image_parameters[ms]['cellsize_in_arcsec'] = [floor(pixscale*1000 + 0.5)/1000] * 2  # in arcsec
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
                field_image_parameters[band][field]['phasecenter'].append(vis_image_parameters[cvis]['phasecenter'])

    return {'visBased':vis_image_parameters, 'fieldBased':field_image_parameters,'filesPerField': files_per_field}
