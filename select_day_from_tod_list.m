function[value]=select_day_from_tod_list(tod_list,day,mydates)
if ~exist('mydates')
  mydates=get_tod_date_from_name (tod_list);
end
fwee=mat2cell (mydates,ones(size(mydates,1),1),size(mydates,2));
ind=strcmp(fwee,day);
value=tod_list(ind);
