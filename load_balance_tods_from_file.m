function[tod_times,myids,node_times]=load_balance_tods_from_file(tod_names,fname,nproc)
tod_times=zeros(size(tod_names));

[names_in,times_in]=read_tod_times(fname);
names_in=get_tod_tags_from_names(names_in);
tags=get_tod_tags_from_names(tod_names);

tt(length(tags),length(tags{1}))=' ';
for j=1:length(tags),
  tt(j,:)=tags{j};
end

nn(length(names_in),length(names_in{1}))=' ';
for j=1:length(names_in),
  nn(j,:)=names_in{j};
end

ntag=length(tags);
for j=1:ntag,
  ind=strmatch(tt(j,:),nn);
  if numel(ind)==1,
    tod_times(j)=times_in(ind);
  end
end

myids=zeros(ntag);
node_times=zeros(nproc,1);


if (1)
  [times_use,ind_use]=sort(tod_times);
  for j=ntag:-1:1,
    ii=ind_use(j);
    assert(times_use(j)==tod_times(ii));
    [a,b]=min(node_times);
    myids(ii)=b;
    node_times(b)=node_times(b)+times_use(j);
  end
else
  myids=zeros(size(tod_names));
  for j=1:ntag,
    [a,b]=min(node_times);
    myids(j)=b;
    node_times(b)=node_times(b)+tod_times(j);
  end
end

mdisp(sprintf('mean node time is %8.3f, min is %8.3f, and max is %8.3f', mean(node_times),min(node_times),max(node_times)));

function[names_in,times_in]=read_tod_times(fname)

lines_in=read_lines(fname);
names_in=cell(size(lines_in));
times_in=zeros(size(lines_in));

nline=length(times_in);
for j=1:nline,
  [a,b]=strtok(lines_in{j});
  [b,c]=strtok(b);
  names_in(j)={a};
  times_in(j)=str2num(b);
end

