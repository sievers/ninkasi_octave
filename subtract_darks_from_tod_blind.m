function[value]=subtract_darks_from_tod_blind(tod)
%subtract off the dark detectors from the tod data.  Don't do anything fancy - just take all the dark 
%detector in a row and subtract off the average of 'em.

data=get_tod_data(tod);
if isempty(data)
  warning('tod is empty in subtract_darks_from_tod_blind');
  return;
end

[row,col]=get_tod_rowcol(tod);col=col+1; %make us unit-offset

darkdat=get_darkdet_avg_cols (tod);
if size(darkdat,2)<max(col),
  darkdat(1,max(col))=0;
end

data=data-darkdat(:,col);%subtract the darks in my column
push_tod_data(data,tod);



