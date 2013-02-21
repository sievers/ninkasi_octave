function[map,ringnest]=read_healpix_map(fname)
fid=fopen(fname,'r');
[aa,bb]=read_fits_header (fid);
[cc,dd]=read_fits_header(fid);

crap=read_fits_bintable (cc,fid);
for j=1:length(cc)
  if strcmp(cc{j,1},'ORDERING')
    ringnest=cc{j,2};
  end
end
map=crap.TEMPERATURE;
fclose(fid);