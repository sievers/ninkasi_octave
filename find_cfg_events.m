function[inds]=find_cfg_events(events,tag)

ii=false(size(events));
nn=length(tag);
for j=1:length(events),
  if strncmp(events{j}.name,tag,nn)
    ii(j)=true;
  end
end

inds=find(ii);