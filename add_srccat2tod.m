function[value]=add_srccat2tod(tod,srccat)
nsrc=numel(srccat.ra);
if isfield(srccat,'oversamp')
  oversamp=srccat.oversamp;
else
  oversamp=0;
end

for j=1:nsrc,
  add_src2tod(tod,srccat.ra(j),srccat.dec(j),srccat.amps(j),srccat.beam,srccat.dx*pi/180,oversamp);
end
