function[map]=set_skymap_from_file_cea(map,fname)

assert(class(map)=='int64');
[fmap,header,names]=fits_image_read(fname,true);
radelt=get_keyword('CDELT1',names,header);
rapix=get_keyword('CRPIX1',names,header);
decdelt=get_keyword('CDELT2',names,header);
decpix=get_keyword('CRPIX2',names,header);
pv=get_keyword('PV2_1',names,header);

nra=get_keyword('NAXIS1',names,header);
ndec=get_keyword('NAXIS2',names,header);



disp([rapix radelt decpix decdelt pv nra ndec]);


set_skymap_cea_predef_c(map,radelt,decdelt,rapix,decpix,pv,ndec,nra);