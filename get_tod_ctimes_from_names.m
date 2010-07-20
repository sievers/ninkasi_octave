function[tags]=get_tod_ctimes_from_names(tod_names)
%return the extension part of a set of TOD names

if ischar(tod_names)
  tod_names={tod_names};
end

if length(tod_names)==1
  tags=parse_tod_name(tod_names{1});
else
  tags=zeros(length(tod_names),1);
  for j=1:length(tod_names)
    tags(j)=parse_tod_name(tod_names{j});
  end
end


function[tag]=parse_tod_name(tod_name)
if tod_name(end)=='/'
  tod_name=tod_name(1:end-1);
end
ind=max(find(tod_name=='/'));
tag=tod_name(ind+1:end);
fwee=min(find(tag=='.'));
tag=str2num(tag(1:fwee-1));