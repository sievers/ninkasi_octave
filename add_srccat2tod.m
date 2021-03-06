function[value]=add_srccat2tod(tod,srccat,varargin)
if iscell(srccat)
  for j=1:numel(srccat),
    add_srccat2tod(tod,srccat{j},varargin{:});
  end
  return
end

nsrc=numel(srccat.ra);
if isfield(srccat,'oversamp')
  oversamp=srccat.oversamp;
else
  oversamp=0;
end
if (1)
  add_srcvec2tod(tod,srccat.ra,srccat.dec,srccat.amps,srccat.beam,srccat.dx*pi/180,oversamp);
else
  for j=1:nsrc,
    %mdisp(['src stuff is: ' num2str([srccat.ra(j) srccat.dec(j) srccat.amps(j) oversamp srccat.dx])]);
    add_src2tod(tod,srccat.ra(j),srccat.dec(j),srccat.amps(j),srccat.beam,srccat.dx*pi/180,oversamp);
  end
end
