function[bands,edges]=get_band_edges_auto(tod,myopts)

%scale_facs=get_struct_mem(myopts,'noise_scale_facs',1);
frac=get_struct_mem(myopts,'bin_frac',0.05); %default to 5% wide bins in noise
min_width=get_struct_mem(myopts,'min_bin_width',3); %default to at least 3 elements per bin
ndata=get_tod_ndata(tod);
dt=get_tod_dt(tod);
edges=zeros(ndata,1);
icur=1;
edges(1)=1;
while edges(icur)<ceil( ndata/2-1)
  delt=ceil((frac)*edges(icur));
  if delt<min_width
    delt=min_width;
  end
  icur=icur+1;
  edges(icur)=edges(icur-1)+delt;
end
edges=edges(1:icur);
if iseven(ndata)
  nn=(ndata/2)+1;
else
  nn=(ndata+1)/2;
end
if edges(end)>nn,
  edges(end)=nn;
end
if edges(end)<edges(end-1)+min_width
  edges=edges([1:end-2 end]);
end

bands=edges/(get_tod_ndata(tod)*get_tod_dt(tod));
