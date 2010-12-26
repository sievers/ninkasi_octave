function[value]=fit_powlaw_plus_offset(data, varargin)

slope_guess=get_keyval_default('slope',-1.5,varargin{:});



datft=fft_r2c(data);
datft=abs(datft(2:end-1));  %cut out the first (since means are undefined) and last (stats depend on even/odd) samples
x=(1:length(datft))';
x=x/x(end);

vecs=[ones(size(x)) x.^slope_guess];
p=linfit(vecft,vecs);

