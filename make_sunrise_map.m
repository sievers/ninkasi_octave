function[offs,sunrise_use,tod_ctimes]=make_sunrise_map(tods,map,fname)

if ~exist('fname')
  fname='/project/rbond/sievers/act/miscellaneous/sunrises_2008.txt';
end


fid=fopen(fname,'r');
sunrise=fscanf(fid,'%d',[1 inf]);
fclose(fid);

tod_ctimes=zeros(size(tods));
sunrise_use=zeros(size(tods));
for j=1:length(tods),
  tod_ctimes(j)=get_tod_ctime_c(tods(j));
  [a,b]=min(abs(tod_ctimes(j)-sunrise));
  sunrise_use(j)=sunrise(b);
end
offs=tod_ctimes-sunrise_use;


if max(abs(offs))>86400/2
  warning('I think I missed some sunrises in make_sunrise_map');
end



make_timemap_octave(tods,map,sunrise_use);
