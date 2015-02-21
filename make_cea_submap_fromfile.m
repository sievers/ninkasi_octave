function[submap,params_new]=make_cea_submap_fromfile(fname,ramin,ramax,decmin,decmax,do_copy)
[params,mm]=get_fits_projection_params(fname);
ra=[ramin ramin ramax ramax];
dec=[decmin decmax decmin decmax];
[xpix,ypix]=radec2pix_fits(ra,dec,params)


xmin=min(floor(xpix))-1;
xmax=max(ceil(xpix))+1;

ymin=min(floor(ypix))-1;
ymax=max(floor(ypix))+1;
params_new=params;
params_new.raoff=params.raoff-xmin;
params_new.decoff=params.decoff-ymin;
mynx=xmax-xmin+1;
myny=ymax-ymin+1;


if ((xmin>=0)&&(xmax<size(mm,1))&&(ymin>=0)&(ymax<size(mm,2)))
  submap=mm((xmin:xmax)+1,(ymin:ymax)-1);
else
  submap=zeros(xmax-xmin+1,ymax-ymin+1);
  sz=size(mm);
  [x1_big,x1_small]=inbound_lims(xmin+1,1,sz(1));
  [y1_big,y1_small]=inbound_lims(ymin+1,1,sz(2));
  [x2_big,x2_small]=inbound_lims(xmax+1,xmax-xmin+1,sz(1));
  [y2_big,y2_small]=inbound_lims(ymax+1,ymax-ymin+1,sz(2));
  submap(x1_small:x2_small,y1_small:y2_small)=mm(x1_big:x2_big,y1_big:y2_big);
end





function[big,small]=inbound_lims(big,small,xmax)
if big<1
  dx=1-big;
  big=1;
  small=1+dx;
end
if big>=xmax
  dx=big-xmax;
  big=xmax;
  small=small-dx;
end
return


