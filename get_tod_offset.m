function[dx,dy]=get_tod_offset(tod_name,offsets)
if iscell(tod_name),
  dx=zeros(length(tod_name),1);
  dy=zeros(length(tod_name),1);
  for j=1:length(tod_name),
    [ddx,ddy]=get_tod_offset(tod_name{j},offsets);
    dx(j)=ddx;
    dy(j)=ddy;
  end
  return
end

crap=strsplit(tod_name,'/',true);
crud=strsplit(crap{end},'.',true);
myct=str2num(crud{1});
ii=find(myct==offsets.ctime);

%set offsets to nan if ctime isn't found
if isempty(ii)
  dx=nan;
  dy=nan;
  return
end
if numel(ii)~=1
  sprintf('extreme oddness in get_tod_offset on %d.',myct)
  ii=ii(1);
end

assert(numel(ii)==1);
dx=offsets.dx(ii);
dy=offsets.dy(ii);
return
