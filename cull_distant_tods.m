function[tods]=cull_distant_tods(tods,opts,dofree)
%cut tods too distant from a given point.  the cut positions may be 
%vectorized, in case one wishes to cover a region.  if cutthresh may be
%either a vector or a scalar.  if a vector, must match the length of RA, and
%each ra will have a separate threshold checked.

if (~isfield(opts,'distthresh'))|(~isfield(opts,'ra'))|(~isfield(opts,'dec'))
  warning('problem in cull_distant_tods.  fields missing, so skipping.');
end

if ~exist('dofree')
  dofree=false;
end


ra=opts.ra;
dec=opts.dec;
thresh=opts.distthresh;

nra=length(ra);
ndec=length(dec);
nthresh=length(thresh);
assert(nra==ndec);
assert( (nthresh==nra)|(nthresh==1));
if nthresh==1,
  thresh=repmat(thresh,size(nra));
end


ind=true(size(tods));

for j=1:nra,
  dists=get_tod_dist(ra(j),dec(j),tods);
  ind(dists>thresh(j))=false;
end

if dofree
  for j=1:length(tods),
    if ~ind(j),
      erase_tod_c(tods(j));          
    end
  end
end
tods=tods(ind);

