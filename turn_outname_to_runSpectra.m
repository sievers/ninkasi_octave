function[fname]=turn_outname_to_runSpectra(fname)
if sum(fname=='/')
  fname=fname(max(find(fname=='/'))+1:end);
end
target='myset_';
ind=findstr(fname,target);
ll=length(target);
fname=[fname(1:ind+ll-1) 'DIVNUM' fname(ind+ll+1:end)];
fname=[fname 'ITERNUM.fits'];
