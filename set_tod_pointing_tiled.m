function[value]=set_tod_pointing_tiled(tods,varargin)
delta=get_keyval_default('delta',0.002,varargin{:});
time_order=get_keyval_default('order',2,varargin{:});
ra_wrap=get_keyval_default('ra_wrap',pi,varargin{:});
for j=1:length(tods),
  [a,b,c,d,e,f]=get_patchy_pointing_fit(tods(j),delta,time_order,ra_wrap);
  set_tod_pointing_tiled_c(tods(j),a,b,c,d,fliplr(e),fliplr(f));
end
