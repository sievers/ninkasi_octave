function[tod_names,use_list]=remove_tods_from_blacklist(tod_names,blacklist)
if ~iscell(blacklist)
  blacklist=read_text_file_comments(blacklist);
end

tod_names_small=get_ctime_part(tod_names);
blacklist=get_ctime_part(blacklist);

use_list=true(size(tod_names));
n1=length(tod_names);
n2=length(blacklist);
for j=1:n1,
  for k=1:n2,
    if (strcmp(tod_names_small{j},blacklist{k}))
      use_list(j)=false;
    end
  end
end

tod_names=tod_names(use_list);


function[names]=get_ctime_part(names)
for j=1:length(names),
  nm=names{j};
  if nm(end)=='/',
    nm=nm(1:end-1);
  end

  ind=1:length(nm);
  ii=ind(nm=='/');
  if length(ii>0),
    nm=nm(max(ii)+1:end);
  end
  names(j)={nm};
end
