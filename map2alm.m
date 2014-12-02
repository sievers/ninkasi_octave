function[alms]=map2alm(map,wt,lmax,mmax)
nside=sqrt(length(map)/12);
assert(nside==round(nside));
if ~exist('wt')
  wt=ones(size(map));
end
if ~exist('lmax')
  lmax=3*nside;
end
if ~exist('mmax')
  mmax=lmax;
end
if min(size(map))==1,
  alms=map2alm_intensity_c(map,lmax,mmax,wt);
else
  error(['it looks like you might have sent in a polarization map with size ' num2str(size(map)) '.  That is not yet supported in map2alm.']);
end

