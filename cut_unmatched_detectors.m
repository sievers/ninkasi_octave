function[asdf]=cut_unmatched_detectors(tods,varargin)
for j=1:length(tods),
  tod=tods(j);


  [rr,cc]=get_tod_rowcol(tod);
  
  match=match_rowcol_pairs(rr,cc,varargin{:});
  for j=1:length(match),
    if match(j)<0,
      cut_detector_c(tod,rr(j),cc(j));
    end
  end
  purge_cut_detectors(tod);
end






