function[newmap]=make_cea_submap(map,ramin,ramax,decmin,decmax,do_copy)
if ~exist('do_copy')
  do_copy=false;
end


%[rapix,decpix,radelt,decdelt,pv,nx,ny]=get_skymap_cea_params_c(map);
[rapix,decpix,radelt,decdelt,pv]=get_skymap_cea_params_c(map);
params.radelt=radelt;
params.decdelt=decdelt;
params.raoff=rapix;
params.decoff=decpix;
params.pv=pv;
params.fittype='cea';

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

[mynx myny]

newmap=allocate_ninkasi_skymap(pi/180,ramin-1e-2,ramax+1e-2,decmin-1e-2,decmax+1e-2);

set_skymap_cea_predef_c(newmap,params_new.radelt,params_new.decdelt,params_new.raoff,params_new.decoff,params_new.pv,mynx,myny);

if (do_copy)
  mm=skymap2octave(map);
  if (1)
    if ((xmin>=0)&&(xmax<size(mm,1))&&(ymin>=0)&(ymax<size(mm,2)))
      submap=mm((xmin:xmax)+1,(ymin:ymax)+1);
    else      
      submap=zeros(xmax-xmin+1,ymax-ymin+1);
      sz=size(mm);
      [x1_big,x1_small]=inbound_lims(xmin+1,1,sz(1));
      [y1_big,y1_small]=inbound_lims(ymin+1,1,sz(2));
      [x2_big,x2_small]=inbound_lims(xmax+1,xmax-xmin+1,sz(1));
      [y2_big,y2_small]=inbound_lims(ymax+1,ymax-ymin+1,sz(2));
      submap(x1_small:x2_small,y1_small:y2_small)=mm(x1_big:x2_big,y1_big:y2_big);
    end
  else
    submap=mm((xmin:xmax)+1,(ymin:ymax)+1);
  end
  octave2skymap(submap,newmap);
end
return


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


