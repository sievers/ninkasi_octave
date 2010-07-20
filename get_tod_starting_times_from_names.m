function[times]=get_tod_starting_times_from_names(tod_names)
%return which second of the day a tod started
times(length(tod_names),1)=0;
for j=1:length(tod_names),
  tt=tod_names{j};
  assert(sum(tt=='.')==2);
  fwee=strtok(tt,'.');
  nn=length(fwee);
  while (nn>1) & (fwee(nn-1)~='/') 
    nn=nn-1;
  end


  if j==1
    t0=str2num(fwee(nn:end));
    fwee=date_from_ctime (str2num(fwee(nn:end)),'HH MM SS');
    [hh,fwee]=strtok(fwee);
    [mm,fwee]=strtok(fwee);
    [ss,fwee]=strtok(fwee);
    assert(isempty(fwee));
    times(j)=str2num(hh)*3600+str2num(mm)*60+str2num(ss);
  else
    tt=str2num(fwee(nn:end));
    dt=tt-t0+times(1);
    times(j)=rem(dt,86400);
    if (times(j)<0)
      times(j)=times(j)+86400;
    end
  end

  %times(j)=str2num(fwee(nn:end));
end

