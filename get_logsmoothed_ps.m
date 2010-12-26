function[smooth_ps]=get_logsmoothed_ps(data,varargin)

smooth_frac=get_keyval_default('frac',0.5,varargin{:});


dataft=fft(data);
nn=floor((size(data,1)+1)/2);
dd=abs(dataft.^2);
zero_val=dd(1);
dd=dd(2:nn,:);
nuvec=(1:size(dd,1))';
lognu=log(nuvec);
dnu=smooth_frac/100;
xvec=lognu(1):dnu:lognu(end);
yvec=interp1(lognu,dd,xvec');

yy=[log(yvec);log(flipud(yvec))];
yvec_smooth=smooth_data_fft(yy,log(smooth_frac/dnu));
yvec_smooth=real(yvec_smooth(1:length(yvec)));

smooth_ps=[zero_val;exp(interp1(xvec,yvec_smooth,lognu,'spline','extrap'))];

%yvec_sm