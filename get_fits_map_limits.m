function[lims]=get_fits_map_limits(mapname,pad)
[params,map,header,names]=get_fits_projection_params(mapname);

if ~exist('pad')
  pad=0;
end

[ra1,dec1]=pix2radec_fits(1-pad,1-pad,params);
[ra2,dec2]=pix2radec_fits(size(map,1)+pad,size(map,2)+pad,params);
ra=[ra1 ra2];
dec=[dec1 dec2];
lims=[min(ra) max(ra) min(dec) max(dec)];