function[tot]=get_observing_time(tods)
tot=0;
for j=1:length(tods),
  tvec=get_tod_tvec(tods(j));
  tot=tot+tvec(end)-tvec(1);
end