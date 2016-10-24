function[tod_names,use_list]=remove_tods_from_blacklist_cbass(tod_names,blacklist)
if ~iscell(blacklist)
  blacklist=read_text_file_comments(blacklist);
end
assert(iscell(blacklist));
if size(blacklist,1)==1,
  blacklist=blacklist';
end

tod_names_small=_get_tod_name_ends(tod_names);

blacklist=_get_tod_name_ends(blacklist);

use_list=true(size(tod_names));
n1=length(tod_names);
n2=length(blacklist);
try
  [all_names,inds]=sort([tod_names_small;blacklist]);
catch
  [all_names,inds]=sort([tod_names_small';blacklist]);
end


for j=1:length(all_names)-1,
  if strcmp(all_names{j},all_names{j+1})
    if inds(j)<n1
      use_list(inds(j))=false;
    end
    if inds(j+1)<n1
      use_list(inds(j+1))=false;
    end
  end
end




%for j=1:n1,
%  for k=1:n2,
%    if (strcmp(tod_names_small{j},blacklist{k}))
%      use_list(j)=false;
%    end
%  end
%end

tod_names=tod_names(use_list);



function[tod_names]=_get_tod_name_ends(tod_names)
for j=1:length(tod_names),
  nm=tod_names{j};
  ii=find(nm=='/');
  if ~isempty(ii)
    nm=nm(ii(end)+1:end);
    tod_names(j)=nm;
  end
end

  



