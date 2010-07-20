function[value]=highpass_map(fname,varargin)
[map,names,vals]=read_map_highpassed(fname,varargin{:});
ind=max(find(fname=='.'));
outname=[fname(1:ind-1)  '_highpassed.fits'];
write_fits(map,outname,vals,names);


