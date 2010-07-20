function[fwee,dates]=get_unique_days_from_tod_list(tod_list)
dates=get_tod_date_from_name (tod_list);
dd=unique(dates,'rows');
fwee=mat2cell (dd,ones(size(dd,1),1),size(dd,2));

