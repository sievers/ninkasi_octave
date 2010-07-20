function[did_i_hit]=does_tod_hit_pixel(tods,map,x,y)
did_i_hit=false(size(tods));
for j=1:length(tods),
  make_weightmap_octave(tods(j),map);
  mm=skymap2octave(map);
  if mm(x,y)>0
    did_i_hit(j)=true;
  end
end

  