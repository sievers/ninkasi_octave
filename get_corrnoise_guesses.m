function[x]=get_corrnoise_guesses(x,tods)
for j=1:length(x.corrnoise),
  crap=get_corrnoise_guess(tods(j).data-skymap2tod(tods(j).ind, x.skymap.map));
  x.corrnoise(j).map=crap.map;
  x.corrnoise(j).vecs=crap.vecs;
end

  