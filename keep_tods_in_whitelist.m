function[tod_names,use_list]=keep_tods_in_whitelist(tod_names,blacklist)
if ~iscell(blacklist)
  blacklist=read_text_file(blacklist);
end


if (1)
  tod_names_small=get_ctime_part(tod_names);
  blacklist=get_ctime_part(blacklist);
  tod_names_small=mycell2complex(tod_names_small);
  blacklist=mycell2complex(blacklist);
else
  %this works, but even though cleaner, is 50% slower
  tod_names_small=names2complex(tod_names);
  blacklist=names2complex(blacklist);
end

use_list=false(size(tod_names));
n1=length(tod_names);
n2=length(blacklist);


for j=1:n1,
  if max(tod_names_small(j)==blacklist)>0
    use_list(j)=true;
  end


  %for k=1:n2,
  %  if (strcmp(tod_names_small{j},blacklist{k}))
  %    use_list(j)=true;
  %  end
  %end
end

tod_names=tod_names(use_list);

function[vals]=names2complex(names)
vals=zeros(size(names));
for j=length(names):-1:1,
  nm=names{j};
  if nm(end)=='/',
    nm=nm(1:end-1);
  end
  crap=strsplit(strtrim(nm),'/');
  crud=strsplit(strtrim(crap{end}),'.');
  vals(j)=str2num(crud{1})+i*str2num(crud{2});
end
  

function[vals]=mycell2complex(names)
vals=zeros(size(names));
for j=1:length(names),
  crud=strsplit(strtrim(names{j}),'.');
  vals(j)=str2num(crud{1})+i*str2num(crud{2});
end

return

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
