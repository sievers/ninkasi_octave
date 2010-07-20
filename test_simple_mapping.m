dt=1/400;
decimate=2;

if exist('decimate')
    dt=dt*decimate;
end

if ~exist('tods')
    [ra1,dec1]=read_pointing('../ninkasi/src/tod_pointing_0_0.dat');
    if exist('decimate'),
        ra1=ra1(1:decimate:end,:);
        dec1=dec1(1:decimate:end,:);
    end
    fac=get_smallest_factor_len(size(ra1,1));
    ra1=ra1(1:fac,:);
    dec1=dec1(1:fac,:);

    [ra2,dec2]=read_pointing('../ninkasi/src/tod_pointing_1_0.dat');
    if exist('decimate'),
        ra2=ra2(1:decimate:end,:);
        dec2=dec2(1:decimate:end,:);
    end
    fac=get_smallest_factor_len(size(ra2,1));
    ra2=ra2(1:fac,:);
    dec2=dec2(1:fac,:);


    minra=min([min(min(ra1)) min(min(ra2))]);
    maxra=max([max(ra1) max(ra2)]);
    mindec=min([min(min(dec1)) min(min(dec2))]);
    maxdec=max([max(dec1) max(dec2)]);
    disp(['limits are ' num2str([minra maxra mindec maxdec]*180/pi)]);
    pix_size=30;
    ra_ind=1+ceil((ra1-minra)*180/pi*3600/pix_size);
    clear ra1;
    dec_ind=1+ceil((dec1-mindec)*180/pi*3600/pix_size);
    clear dec1;

    ra_ind2=1+ceil((ra2-minra)*180/pi*3600/pix_size);
    clear ra2;
    dec_ind2=1+ceil((dec2-mindec)*180/pi*3600/pix_size);
    clear dec2;

    map_nx=max([max(max(ra_ind)) max(max(ra_ind2))]);
    map_ny=max([max(max(dec_ind)) max(max(dec_ind2))]);

    ind1=int32(dec_ind-1+map_ny*(ra_ind-1));
    %ind1=[ra_ind dec_ind];
    clear ra_ind;
    clear dec_ind;

    ind2=int32(dec_ind2-1+map_ny*(ra_ind2-1));
    %ind2=[ra_ind2 dec_ind2];
    clear dec_ind2;
    clear ra_ind2;
    tods(1).ind=ind1;
    clear ind1;
    tods(2).ind=ind2;
    clear ind2;
    tods(1).ndata=size(tods(1).ind,1);
    tods(2).ndata=size(tods(2).ind,1);
    tods(1).ndet=size(tods(1).ind,2);
    tods(2).ndet=size(tods(2).ind,2);

end
tic
weight_map=tod2map(tods(1).ind,map_nx,map_ny)+tod2map(tods(2).ind,map_nx,map_ny);
toc

if (1)
    map=10*generate_cmb_map([map_ny map_nx],pix_size/60);
else
    map=zeros(map_nx,map_ny);
    r1=(1:map_nx)-map_nx/2;
    r2=(1:map_ny)'-map_ny/2;
    rr=repmat(r1.^2,[map_ny 1])+repmat(r2.^2,[1 map_nx]);
    map=exp(-0.5*rr/100);
end
mapvec=reshape(map,[numel(map) 1]);
ntod=2;


tic
for tt=1:ntod,
    data=zeros(size(tods(tt).ind));
    for j=1:size(data,2),
        data(:,j)=mapvec(tods(tt).ind(:,j));
    end
    tods(tt).data_org=data;
    %tods(tt).data=zeros(size(tods(tt).ind));
    %tods(tt).data(:,j)=mapvec(tods(tt).ind(:,j));
end
toc
dirty_map=0;
dirty_map2=0;
for j=1:ntod,
    tods(j).data=tods(j).data_org+0.5*randn(size(tods(j).data_org));
    dirty_map=dirty_map+tod2map(tods(j).data,tods(j).ind,map_nx,map_ny);
    [tods(j).cm,tods(j).skyspec]=get_1overf_noise(size(tods(j).data,1),1,0,4,-1.5,dt,false);
    tods(j).data=tods(j).data+repmat(tods(j).cm,[1 size(tods(j).data,2)]);
    tods(j).cm_guess=mean(tods(j).data,2);
    dirty_map2=dirty_map2+tod2map(tods(j).data-repmat(tods(j).cm_guess,[1 tods(j).ndet]),tods(j).ind,map_nx,map_ny);
end
mm=divide_maps(dirty_map2,weight_map);clf;imagesc(mm);axis equal;colorbar


%dirty_map2=tod2map(data1_noisy,ind1,map_nx,map_ny)+tod2map(data2_noisy,ind2,map_nx,map_ny);

skymap.map=double(map);
for j=1:length(tods),
    corrnoise(j).map=0*tods(j).cm;
    corrnoise(j).vecs=ones(1,tods(j).ndet);
end
mapset.skymap=skymap;
mapset.corrnoise=corrnoise;

clear x;
x.skymap=mapset.skymap;
x.skymap.map=0*x.skymap.map;
%x=mapset;
%x.skymap.map=0*x.skymap.map;
for j=1:length(tods),
    crud=get_corrnoise_guess(tods(j).data);
    crud.prior={tods(j).skyspec};
    crud.has_prior=true;
    x.corrnoise(j)=crud;
end


%code to do super-iteration
%x3=x2;x3.corrnoise(1)=get_corrnoise_guess(tods(1).data-skymap2tod(tods(1).ind,x2.skymap.map));x3.corrnoise(2)=get_corrnoise_guess(tods(2).data-skymap2tod(tods(2).ind,x2.skymap.map));
%precon.skymap.map=weight_map;



