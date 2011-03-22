function[nbolo,nsamp,data,flags,ra,dec,el]=read_raw_bolocam(fname)
fid=fopen(fname,'r');
nbolo=fread(fid,1,'int');
nsamp=fread(fid,1,'int');
data=fread(fid,[nbolo nsamp ],'double')';
flags=fread(fid,[nbolo nsamp],'char')';
ra=fread(fid,[nbolo nsamp ],'double')';
dec=fread(fid,[nbolo nsamp ],'double')';
el=fread(fid,[nsamp 1],'double');
fclose(fid);


8*(numel(ra)+numel(dec)+numel(el)+numel(data))+1*numel(flags)+8



