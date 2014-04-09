function[asd]=restore_input_map()

initialize_mapping_scripts;

fnames=__get_files('ls -1 *input_used.fits | grep -v single');

tails={'200','500','1000'};

nside=512;
pixsize=0.01;
map=allocate_ninkasi_skymap(pixsize,0.1,0.2,0.1,0.2);

set_skymap_healpix_ring_c(map,nside);



for j=1:length(fnames),
  froot=fnames{j};
  mm=read_cbass_intensity_map(froot);
  nn=sqrt(numel(mm)/12);
  if nn~=nside
    disp(['changing nside to ' num2str(nn)]);
    set_skymap_healpix_ring_c(map,nn);
    nside=nn;
  end

  tag='_input_used';
  ii=strfind(froot,tag);
  froot=froot(1:ii-1);
  disp(froot);

  for k=1:length(tails),
    inname=[froot tails{k} '.fits'];
    outname=[froot tails{k} '_restored'];
    if exist(inname)
      disp(['restoring ' inname ' into ' outname]);
      m2=read_cbass_intensity_map(inname);
      m2=m2+mm;
      octave2skymap(m2,map);
      write_map(map,outname);
    end
  end
end
destroy_map(map);

return

m1=read_cbass_intensity_map('../6month_map_512.fits');

froots={
'3month_ml_mapsub_fine_nside_512_500'
};

pixsize=0.01;
map=allocate_ninkasi_skymap(pixsize,0.1,0.2,0.1,0.2);
nside=sqrt(numel(m1)/12);
set_skymap_healpix_ring_c(map,nside);

for j=1:length(froots)
  froot=froots{j};
  
  m2=read_cbass_intensity_map([froot '.fits']);
  
  mm=m1+m2;
  octave2skymap(mm,map);
  write_map(map,[froot '_restored.fits']);
end


function[fnames]=__get_files(ff)
[a,crud]=system(ff);
nn=sprintf('\n');
fnames=strsplit(strtrim(crud),nn);