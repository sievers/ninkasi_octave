function[hits]=tod_hits_srccat(tod,srccat)
if iscell(srccat),
  hits=cell(size(srccat));
  for j=1:numel(srccat),
    hits(j)={tod_hits_srccat(tod,srccat{j})};
  end
  return
end
nsrc=numel(srccat.ra);
dist=srccat.dx*pi/180*numel(srccat.beam);
hits=zeros(nsrc,1);
for j=1:nsrc,
  hits(j)=tod_hits_source_c(srccat.ra(j),srccat.dec(j),dist,tod);
end


