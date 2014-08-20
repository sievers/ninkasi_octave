function[map]=expand_periodic_map(map,targ)
sz_targ=size(targ);
sz=size(map);
assert(numel(sz)==numel(sz_targ));
if numel(sz)==3,
  while size(map,2)<size(targ,2)
    map=repmat(map,[1 2 1]);
  end
  map=map(:,1:size(targ,2),:);
  while size(map,3)<size(targ,3)
    map=repmat(map,[1 1 2]);
  end
  map=map(:,:,1:size(targ,3));
else
  assert(numel(sz)==2)
  while size(map,1)<size(targ,1)
    map=repmat(map,[ 2 1]);
  end
  map=map(1:size(targ,2),:);
  while size(map,2)<size(targ,2)
    map=repmat(map,[ 1 2]);
  end
  map=map(:,1:size(targ,2));
end

