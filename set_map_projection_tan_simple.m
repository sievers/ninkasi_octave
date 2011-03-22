function[value]=set_map_projection_tan_simple(map,varargin)
[ramin,ramax,decmin,decmax,pixsize]=get_skymap_rect_params_c (map);



racent=0.5*(ramin+ramax);
racent=get_keyval_default('racent',racent,varargin{:});
deccent=0.5*(decmin+decmax);
deccent=get_keyval_default('deccent',deccent,varargin{:});
pixsize=get_keyval_default('pixsize',pixsize,varargin{:});


rapix_approx=0.5*(ramax-ramin)*cos(deccent)/pixsize;
decpix_approx=0.5*(decmax-decmin)/pixsize;

mdisp(['ra/dec cents are ' num2str([racent deccent])])

mdisp([pixsize rapix_approx decpix_approx racent deccent])
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

[rax' decx']
xmin=min(rax);
xmax=max(rax);
ymin=min(decx);
ymax=max(decx);
pad=2;

nra=ceil(xmax-xmin)+2*pad;
ndec=ceil(ymax-ymin)+2*pad;
rapix=rapix_approx-xmin+pad;
decpix=decpix_approx-ymin+pad;
[pixsize rapix decpix racent deccent nra ndec]
set_skymap_tan_predef_c (map,pixsize,rapix,decpix,racent,deccent,nra,ndec);





