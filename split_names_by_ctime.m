function[names_low,names_high]=split_names_by_ctime(tod_names,round_up)
ctimes=get_tod_ctimes_from_names(tod_names);

names_low=tod_names(ctimes<median(ctimes));
names_high=tod_names(ctimes>median(ctimes));
if ~iseven(length(tod_names)),
  if round_up,
    names_high=[names_high tod_names(ctimes==median(ctimes))];
  else
    names_low=[names_low tod_names(ctimes==median(ctimes))];
  end
end
