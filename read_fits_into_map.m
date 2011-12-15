function[mapptr,params]=read_fits_into_map(fname)
[params,map]=get_fits_projection_params(fname);

pixsize=30*pi/180/3600;
lims=[0 0.1 0 0.1]*pi/180;
mapptr=allocate_ninkasi_skymap(pixsize,lims(1)-1e-2,lims(2)+1e-2, lims(3)-1e-2,lims(4)+1e-2);
switch(params.fittype)
 case{'cea'}
  set_skymap_cea_predef_c(mapptr,params.radelt,params.decdelt,params.raoff,params.decoff,params.pv,size(map,1),size(map,2));
  octave2skymap(map,mapptr);
 otherwise
  error(['Currently unsupported projection in read_fits_into_map: ' params.fittype]);
end



