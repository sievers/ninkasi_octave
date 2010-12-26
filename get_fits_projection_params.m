function[params,map]=get_fits_projection_params(fname)

[map,header,names]=fits_image_read(fname);
fittype=get_keyword_string('CTYPE1',names,header);
fittype(fittype=='''')=' ';
fittype=strtrim(fittype);
switch lower(fittype)
  case {'ra---cea'}
   params=parse_cea_params(names,header);
   
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
