function[calfacs]=parse_tod_recal_file(fname)
%if pass in a character string, we need to read teh file, otherwise it has already been read
fid=fopen(fname,'r');
vec=fread(fid,inf,'char=>char');fclose(fid);
vec(vec=='a')=' ';
vec(vec=='r')=' ';
%vec(vec=='.')=' ';
vec(vec=='e')=' ';
vec(vec=='+')=' ';

vec(vec==sprintf('\n'))=' ';
ii=find(vec=='.');
vec(ii(1:3:end))=' ';
vec(ii(2:3:end))=' ';


vals=sscanf(vec','%f');

ct1=floor(vals(1:5:end));
ct2=vals(2:5:end);
arr=vals(3:5:end);
cal=vals(4:5:end);
cal=cal.*10.^vals(5:5:end);

calfacs.ct1=ct1;
calfacs.ct2=ct2;
calfacs.arr=arr;
calfacs.cal=cal;
