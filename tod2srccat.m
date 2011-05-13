function[srccat]=tod2srccat(tod,srccat)
nsrc=numel(srccat.ra);
if isfield(srccat,'oversamp')
  oversamp=srccat.oversamp;
else
  oversamp=0;
end
srccat.amps=srccat.amps+tod2srcvec(tod,srccat.ra,srccat.dec,srccat.beam,srccat.dx*pi/180,oversamp);

