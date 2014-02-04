function[mymap]=read_cbass_intensity_map(fname)
fid=fopen(fname,'r');
h=read_fits_header(fid);
h2=read_fits_header(fid);
asdf=read_fits_bintable_test(h2,fid);
mymap=asdf.TEMPERATURE;
mymap=reshape(mymap',[numel(mymap) 1]);
fclose(fid);
