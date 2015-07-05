function[tod_times,myids,node_times]=load_balance_tods_from_file(tod_names,fname,nproc,varargin)
%how to assign TODs with missing time info.
%options include mean, median, max,min,cut
%default of min means the load balancing happens as normal,
%then they will get dealt out at the end.
do_missing=get_keyval_default('missing','min',varargin{:});

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

if sum(tod_times>0)==0,
  warning('did not find any tods with time information.')
  tod_times=ones(size(tod_names));
  myids=(1:numel(tod_names))';
  myids=rem(myids-1,nproc)+1;
  node_times=ones(size(tod_names));
  return
end


if sum(tod_times==0)>0,
  mdisp(['warning - in load_balance_tods_from_file, have TODs with missing time info']);
  switch(do_missing)
   case{'min'}
    tod_times(tod_times==0)=0.9999*min(tod_times(tod_times>0));
   case{'mean'}
    tod_times(tod_times==0)=mean(tod_times(tod_times>0));
   case{'median'} 
    tod_times(tod_times==0)=median(tod_times(tod_times>0));
   case{'max'}
    tod_times(tod_times==0)=max(tod_times(tod_times>0));
   case{'cut'}
    ;  %don't run TODs with no time info
   case{'ignore'}
    ;  
   otherwise
    mdisp(['unrecognized input ' do_missing ' for missing TOD time info in load_balance_tods_from_file.  Assuming cut.']);
  end
end



myids=zeros(ntag,1);
node_times=zeros(nproc,1);


if (1)
  [times_use,ind_use]=sort(tod_times);

  if (min(times_use)==0)&(strcmp(do_missing,'ignore')==0)
    imin=max(find(times_use==0))+1;
  else
    imin=1;
  end
  for j=ntag:-1:imin,
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

