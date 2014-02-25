function[value]=align_cea_maps(maplist,outlist,varargin)
fac=get_keyval_default('scale_fac',1.0,varargin{:});

mpi_init;
nmap=length(maplist);
for j=1:nmap,
  [myproj,mymap]=get_fits_projection_params(maplist{j});
  if (fac~=1)
    mymap=mymap*fac;
  end
  [myra1,mydec1]=pix2radec_fits(1,1,myproj);
  [myra2,mydec2]=pix2radec_fits(size(mymap,1),size(mymap,2),myproj);
  lims(j,:)=[min(myra1,myra2) max(myra1,myra2) min(mydec1,mydec2) max(mydec1,mydec2)];
  proj(j)={myproj};
  maps(j)={mymap};
  if j==1
    pv_targ=myproj.pv;
    pixtarg=myproj.radelt;,
    pixtarg2=myproj.decdelt;,
  else
    assert(myproj.pv==pv_targ);
    assert(myproj.radelt==pixtarg);
    assert(myproj.decdelt==pixtarg2);
  end

end

pixsize=abs(pixtarg)
lims

alllims=[min(lims(:,1)) max(lims(:,2)) min(lims(:,3)) max(lims(:,4))]
alllims2=alllims+pixsize*[-2 2 -2 2];
mapptr=allocate_ninkasi_skymap(abs(pixtarg),alllims2(1)*pi/180,alllims2(2)*pi/180,alllims2(3)*pi/180,alllims2(4)*pi/180);

set_skymap_cea_simple_predef(mapptr,abs(pixtarg),pv_targ);
mm=skymap2octave(mapptr);
write_map(mapptr,outlist{1});
[outproj,mm2]=get_fits_projection_params([outlist{1} '.fits']);

for j=1:nmap,
  mymap=maps{j};
  mm=0*mm;
  if (1)
    [x1,y1]=radec2pix_fits(lims(j,1),lims(j,3),outproj);
    [x2,y2]=radec2pix_fits(lims(j,2),lims(j,4),outproj);
    
    [xx,yy]=radec2pix_fits(alllims(2),alllims(3),outproj);

    xmin=round(min(x1,x2));
    ymin=round(min(y1,y2));
    xmax=xmin+size(mymap,1)-1;
    ymax=ymin+size(mymap,2)-1;
    disp([xmin xmax ymin ymax,xx,yy])
    mm(xmin:(xmin+size(mymap,1)-1),ymin:(ymin+size(mymap,2)-1))=mymap;

  else
    %this code is broken because set_skymap_cea_simple_predef doesn't always seem to do exactly what I think it should
    decoff=round( (lims(j,3)-alllims(3))/pixsize);  
    raoff=-round((lims(j,2)-alllims(2))/pixsize);
    ramax=raoff+size(mymap,1);
    decmax=decoff+size(mymap,2);
    disp([raoff ramax decoff decmax])
    
    mm(raoff+1:ramax,decoff+1:decmax)=mymap;
  end
  octave2skymap(mm,mapptr);
  write_map(mapptr,outlist{j});
end


