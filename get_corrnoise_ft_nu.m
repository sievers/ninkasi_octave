function[vecft,x]=get_corrnoise_ft_nu(corrnoise,varargin)
%Get the power spectrum and x coordinates of corrnoise for use in priors.

do_abs=get_keyval_default('do_abs',true,varargin{:});


if do_abs
  vecft=abs(fft_r2c(corrnoise));
  vecft=vecft(2:end-1,:);
  x=(1:length(vecft))';
  x=x/max(x);
else
  vecft=(fft_r2c(corrnoise));
  x=(0:length(vecft)-1)';
  x=x/(max(x)-1);
end

