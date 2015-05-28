function[tags,tags2]=get_tod_ctimes_from_names(tod_names)
%return the extension part of a set of TOD names

if ischar(tod_names)
  tod_names={tod_names};
end

if length(tod_names)==1
  [tags,tags2]=parse_tod_name(tod_names{1});
else
  tags=zeros(length(tod_names),1);
  tags2=zeros(length(tod_names),1);
  for j=1:length(tod_names)
    [tags(j),tags2(j)]=parse_tod_name(tod_names{j});
    %[tags2(j)]=parse_tod_name(tod_names{j});
  end
end


function[tag,tag2]=parse_tod_name(tod_name)
if tod_name(end)=='/'
  tod_name=tod_name(1:end-1);
end
if (0)
  ind=max(find(tod_name=='/'));
  tag=tod_name(ind+1:end);
  fwee=min(find(tag=='.'));

else
  [aa,tag,bb]=fileparts(tod_name);
  %fwee=min(find(tag=='.'));
  %tag=str2num(tag(1:fwee-1));
  tags=strsplit(strtrim(tag),'.');
  tag=str2num(tags{1});
  tag2=str2num(tags{2});

end
