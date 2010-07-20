function[params]=set_cea_from_file(fname,map,varargin) 
scale_fac=get_keyval_default('scale_fac',1,varargin{:});
[inmap,header,names]=fits_image_read(fname);
if isempty(inmap)
  disp(['unable to read map ' fname ' from disk in set_cea_from_file.']);
  return
end
if scale_fac~=1
  inmap=inmap*scale_fac;
end

radelt=get_keyword('CDELT1',names,header,0);
decdelt=get_keyword('CDELT2',names,header,0);

rapix=get_keyword('CRPIX1',names,header,0);
decpix=get_keyword('CRPIX2',names,header,0);

pv=get_keyword('PV2_1',names,header,0);
nra=get_keyword('NAXIS1',names,header,0);
ndec=get_keyword('NAXIS2',names,header,0);

%printf('args are %14.8f %14.8f %d %d %14.8f  %d %d\n',radelt,decdelt,rapix,decpix,pv,nra,ndec);

if strcmp(class(map),'int64')
  set_skymap_cea_predef_c(map,radelt,decdelt,rapix,decpix,pv,nra,ndec);
  octave2skymap_c(inmap,map);
else
  assert(isstruct(map));
  assert(isfield(map,'mapptr'));
  set_skymap_cea_predef_c(map.mapptr,radelt,decdelt,rapix,decpix,pv,nra,ndec);
  map.skymap=inmap;
  octave2skymap_c(inmap,map.mapptr);
end