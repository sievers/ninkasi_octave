function[cal]=find_tod_recal_fac(tod_name,recal_facs)
%recal facs should be gotten via parse_tod_recal_file.m
%from one of Matthew Hasselfield's recal files

ii=find(tod_name=='/');
if ~isempty(ii)
  tod_name=tod_name(ii(end)+1:end);
end
disp(tod_name)
tod_name=strtrim(tod_name); %just in case there's trailing whitesapce
ii=find(tod_name=='.')
ct1=str2num(tod_name(1:ii(1)-1));
ct2=str2num(tod_name(ii(1)+1:ii(2)-1));
arr=str2num(tod_name(end));
ii=(ct1==recal_facs.ct1)&(ct2==recal_facs.ct2)&(arr==recal_facs.arr);
assert(sum(ii)==1)
cal=recal_facs.cal(ii);

