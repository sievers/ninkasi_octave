function[params,map,header,names]=get_fits_projection_params(fname)

[map,header,names]=fits_image_read(fname);
fittype=get_keyword_string('CTYPE1',names,header);
fittype(fittype=='''')=' ';
fittype(fittype==0)=' ';
fittype=strtrim(fittype);

strcmp(lower(fittype),'ra---cea')
switch lower(fittype)
  case {'ra---cea'}
   params=parse_cea_params(names,header);
 case {'ra---tan'}
  params=parse_tan_params(names,header);
 otherwise
  warning(['unsupported type - ' fittype]);
end








function[params]=parse_cea_params(names,header);
params.fittype='cea';
params.radelt=get_keyword('CDELT1',names,header);
params.raoff=get_keyword('CRPIX1',names,header);
params.decdelt=get_keyword('CDELT2',names,header);
params.decoff=get_keyword('CRPIX2',names,header);
params.pv=get_keyword('PV2_1',names,header);

function[params]=parse_tan_params(names,header);
params.fittype='tan';
params.radelt=get_keyword('CD1_1',names,header);
params.raoff=get_keyword('CRPIX1',names,header);
params.decdelt=get_keyword('CD2_2',names,header);
params.decoff=get_keyword('CRPIX2',names,header);
params.racent=get_keyword('CRVAL1',names,header);
params.deccent=get_keyword('CRVAL2',names,header);
