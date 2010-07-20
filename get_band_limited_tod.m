function[tod]=get_band_limited_tod(tod,nu0,nu1,varargin)
dt=get_keyval_default('dt',1/200,varargin{:});
n=size(tod,1);
nu=(0:n-1)'/(n*dt);
ind=(nu>nu0)&(nu<nu1);
ind(ceil((end+1)/2)+1:end)=0;
ind=ind|[0;flipud(ind(2:end))];

todft=fft(tod);
todft(~ind,:)=0;
tod=ifft(todft);
