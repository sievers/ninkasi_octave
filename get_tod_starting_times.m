function[times]=get_tod_starting_times(tods)
times=zeros(size(tods));
for j=1:length(tods),
  fwee=date_from_ctime (get_tod_ctime_c(tods(j)),'HH MM SS');
  [hh,fwee]=strtok(fwee);
  [mm,fwee]=strtok(fwee);
  [ss,fwee]=strtok(fwee);
  assert(isempty(fwee));
  times(j)=str2num(hh)*3600+str2num(mm)*60+str2num(ss);
end
