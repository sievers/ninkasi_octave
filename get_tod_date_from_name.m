function[dates]=get_tod_date_from_name(tod_names)

dates(length(tod_names),1)=' ';
for j=1:length(tod_names),
  fwee=strtok(tod_names{j},'.');
  ind=1:length(fwee);
  crap=ind(fwee=='/');
  mydate=fwee(crap(end-1)+1:crap(end)-1);
  dates(j,1:length(mydate))=mydate;
end

  