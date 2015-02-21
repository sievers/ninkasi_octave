function[map]=read_map_into_ptr(fname)
[params,mm]=get_fits_projection_params(fname);
[ra1,dec1]=pix2radec_fits(1,1,params);
[ra2,dec2]=pix2radec_fits(size(mm,1),size(mm,2),params);
pixsize=params.decdelt*pi/180;
ra=[ra1 ra2];
dec=[dec1 dec2];
lims=[min(ra) max(ra) min(dec) max(dec)]*pi/180;

map=allocate_ninkasi_skymap(pixsize,lims(1)-0.5*pixsize,lims(2),lims(3),lims(4)+0.5*pixsize);

set_skymap_cea_simple_predef_c(map,pixsize*180/pi,params.pv);

crud=skymap2octave(map);
octave2skymap(mm,map);
