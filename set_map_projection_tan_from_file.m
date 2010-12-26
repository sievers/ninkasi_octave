function[value]=set_map_projection_tan_from_file(map,fname)
[aa,header,names]=fits_image_read(fname);
ra_cent=get_keyword('CRVAL1',names,header)*pi/180;
dec_cent=get_keyword('CRVAL2',names,header)*pi/180;
rapix=get_keyword('CRPIX1',names,header);
decpix=get_keyword('CRPIX2',names,header);
pv=get_keyword('PV2_1',names,header,0);
radelt=get_keyword('CD1_1',names,header)*pi/180;
decdelt=get_keyword('CD2_2',names,header)*pi/180;
nra=size(aa,1);
ndec=size(aa,2);

set_skymap_tan_explicit_c(map,rapix,decpix,radelt,decdelt,pv,ra_cent,dec_cent,nra,ndec);

return



[ramin,ramax,decmin,decmax,pixsize]=get_skymap_rect_params_c (map);

racent=0.5*(ramin+ramax);
deccent=0.5*(decmin+decmax);

rapix_approx=0.5*(ramax-ramin)*cos(deccent)/pixsize;
decpix_approx=0.5*(decmax-decmin)/pixsize;

disp(['ra/dec cents are ' num2str([racent deccent])])

set_skymap_tan_predef_c (map,pixsize,rapix_approx,decpix_approx,racent,deccent,1,1);
ra_cents=[ramin racent ramax];
dec_cents=[decmin deccent decmax];
ii=0;
for j=1:length(ra_cents),
  for k=1:length(dec_cents),
    ii=ii+1;
    [rax(ii),decx(ii)]=get_pix_from_radec_c (map,ra_cents(j),dec_cents(k));
  end
end
xmin=min(rax);
xmax=max(rax);
ymin=min(decx);
ymax=max(decx);
pad=2;

nra=ceil(xmax-xmin)+2*pad;
ndec=ceil(ymax-ymin)+2*pad;
rapix=rapix_approx-xmin+pad;
decpix=decpix_approx-ymin+pad;
set_skymap_tan_predef_c (map,pixsize,rapix,decpix,racent,deccent,nra,ndec);





