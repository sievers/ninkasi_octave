function[fwee]=select_nth_day_from_tod_list(tod_list,which_day)
[dd,dates]=get_unique_days_from_tod_list(tod_list);
if ischar(which_day)
  fwee=select_day_from_tod_list(tod_list,which_day,dates);
else
  assert(which_day>0);
  assert(which_day <=length(dd));
  fwee=select_day_from_tod_list(tod_list,dd{which_day},dates);
end



