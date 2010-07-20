function[offs,sunrise_use,tod_ctimes]=get_time_from_sunrise(tods,fname)
%get the distance from sunrise of a set of TODs.  Can work either on
%TOD names or actual TODs.  Args are (tods,fname) where fname is the 
%filename containing sunrises.

if ~exist('fname')
  fname='/project/rbond/sievers/act/miscellaneous/sunrises_2008.txt';
end


fid=fopen(fname,'r');
sunrise=fscanf(fid,'%d',[1 inf]);
fclose(fid);


if iscell(tods)
  tod_ctimes=get_tod_ctimes_from_names(tods)';
else
  tod_ctimes=zeros(size(tods));
  for j=1:length(tods),
      tod_ctimes(j)=get_tod_ctime_c(tods(j));
  end
end

sunrise_use=zeros(size(tods));

for j=1:length(tods),
  [a,b]=min(abs(tod_ctimes(j)-sunrise));
  sunrise_use(j)=sunrise(b);
end
offs=tod_ctimes-sunrise_use;


if max(abs(offs))>86400/2
  warning('I think I missed some sunrises in make_sunrise_map');
end
