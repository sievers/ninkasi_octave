function[value]=cut_tod_turnarounds(tods,varargin)
if length(tods)>1
  for j=1:length(tods),
    cut_tod_turnarounds(tods(j),varargin{:});
  end
else
  assert(class(tods)=='int64');
  cut_len=get_keyval_default('cut_len',0.5,varargin{:});
  cut_pre=get_keyval_default('cut_pre',cut_len,varargin{:});
  cut_post=get_keyval_default('cut_post',cut_len,varargin{:});
  
  turns=find_tod_turnarounds (tods);
  dt=get_tod_dt(tods);
  for j=1:length(turns),
    cut_tod_global_c(tods,round(turns(j)-cut_pre/dt),round(turns(j)+cut_post/dt));
  end
end


