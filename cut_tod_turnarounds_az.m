function[big_pre,big_post,turns]=cut_tod_turnarounds_az(tods,varargin)
if length(tods)>1
  for j=1:length(tods),
    cut_tod_turnarounds_az(tods(j),varargin{:});
  end
else
  assert(class(tods)=='int64');
  cut_len=get_keyval_default('cut_len',0.5,varargin{:});
  cut_pre=get_keyval_default('cut_pre',cut_len,varargin{:});
  cut_post=get_keyval_default('cut_post',cut_len,varargin{:});
  
  turns=find_tod_turnarounds (tods);
  cut_pre=cut_pre*pi/180;
  cut_post=cut_post*pi/180;
  [alt,az]=get_tod_altaz(tods);
  ndata=get_tod_ndata(tods);
  big_pre=zeros(size(turns));
  big_post=zeros(size(turns));
  for j=1:length(turns),
    j1=turns(j);
    while (j1>1) & (abs(az(j1)-az(turns(j)))<cut_pre)
      j1=j1-1;   
    end
    j2=turns(j);
    while (j2<ndata) & (abs(az(j2)-az(turns(j)))<cut_post)
      j2=j2+1;
    end
    big_pre(j)=j1;
    big_post(j)=j2;
    cut_tod_global_c(tods,j1,j2);
  end
end


