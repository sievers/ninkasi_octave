function[segs]=inds2segs(inds)
%transform a list of indices into line segments.
if length(inds)==1,
  segs=[inds inds];
  return;
end

df=diff(inds);

ind_end=[find(df>1) length(inds)];
ind_beg=[1 1+ind_end(1:end-1)];

seg_beg=inds(ind_beg);
seg_end=inds(ind_end);

if size(seg_beg,1)==1,
  segs=[seg_beg;seg_end];
else
  segs=[seg_beg seg_end];
end
