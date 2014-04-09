function[value]=apply_wtmap_to_map(map,wt)
mm=skymap2octave(map);
ww=skymap2octave(wt);
mm(ww>0)=mm(ww>0)./ww(ww>0);
octave2skymap(mm,map);

