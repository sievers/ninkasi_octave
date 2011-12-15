function[data,ra,dec,mask,timevec,names]=read_sanepic(fname)
fid=fopen(fname,'r','ieee-be');
[keys,vals]=read_fits_header(fid);
[data,tag]=read_fits_hdu(fid);assert(strcmp(tag,'signal'))
[ra,tag]=read_fits_hdu(fid);assert(strcmp(tag,'ra'))
[dec,tag]=read_fits_hdu(fid);assert(strcmp(tag,'dec'))
[mask,tag]=read_fits_hdu(fid);assert(strcmp(tag,'mask'))
[timevec,tag]=read_fits_hdu(fid);assert(strcmp(tag,'time'))
det_names=read_fits_hdu(fid);names=det_names.names; %this fails if we're confused

fclose(fid);

