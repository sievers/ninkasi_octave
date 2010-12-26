function[fitp]=fit_1overf_noise(vecs,varargin)
vecft=fft(vecs);
nn=floor((size(vecft,1)+1)/2);
vecft=vecft(1:nn,:);

minfreq=get_keyval_default('minfreq',0.05,varargin{:});
maxfreq=get_keyval_default('maxfreq',0.9,varargin{:});

nuvec=(0:nn-1)';
data=abs(vecft.^2);
imin=round(minfreq*nn);
imax=round(maxfreq*nn);
nuvec=nuvec(imin:imax);
data=data(imin:imax,:);

whos
loglog(data(:,1));


