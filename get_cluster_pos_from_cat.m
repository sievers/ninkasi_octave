function[ra,dec,name]=get_cluster_pos_from_cat(name,mycat)
if ~exist('mycat')
  mycat=read_cluster_catalog;
end

ok_ind=[];
for j=1:length(mycat),
  if strcmp(name,mycat(j).name)
    ra=mycat(j).ra;
    dec=mycat(j).dec;
    return
  end
  if strncmp(name,mycat(j).name,length(name))
    ok_ind(end+1)=j;
  end
end
if length(ok_ind)==0
  if name(1)=='J'
    [ra,dec,name]=get_cluster_pos_from_cat(['ACT_' name ],mycat);
    return;
  end
  if ~isempty(str2num(name(1))) %likely just got an RA in here
    [ra,dec,name]=get_cluster_pos_from_cat(['ACT_J' name],mycat);
    return;
  end
  %got something I don't even know how to guess at.
  error(['name ' name ' not found in catalog in get_cluster_pos_from_cat.']);
end
if length(ok_ind)==1
  ra=mycat(ok_ind).ra;
  dec=mycat(ok_ind).dec;
  name=mycat(ok_ind).name;
  return
end
if length(ok_ind)>1
  error(['found too many matches in get_cluster_pos_from_cat on name ' name ]);
end
assert(1==0);  %one should never get here