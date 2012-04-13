function[det_noises,vecs]=read_tod_noise_banded_projvec_octave(fname)
fid=fopen(fname,'r');
if (fid==-1)
  error(['Unable to read noises from ' fname]);
  return;
end
ndet=fread(fid,1,'int');
nband=fread(fid,1,'int');
edges=fread(fid,[1+nband 1],'int');
nvecs=fread(fid,[nband 1],'int');
det_noises=fread(fid,[ndet nband],'double');
vecs=cell(nband,1);
for j=1:nband,
  vecs(j)={fread(fid,[ndet nvecs(j)],'double')};
end
fclose(fid);