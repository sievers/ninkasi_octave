function[fitp]=get_mustang_noise(tod,varargin)

global nuvec
global datasqr



data=get_tod_data_cleaned(tod,varargin{:});
coeffs=zeros(3,size(data,2));
nuvec=get_tod_nuvec(tod);

minfreq=get_keyval_default('minfreq',0.5*nuvec(1),varargin{:});
maxfreq=get_keyval_default('maxfreq',nuvec(end),varargin{:});

dataft=fft_r2c(data);
nuvec=nuvec(1:size(dataft,1));
ind=(nuvec>minfreq)&(nuvec<maxfreq);
ind(end)=false;  %don't use the last sample, since much care must be taken if tods are even/odd length
clear data;

dataft=dataft(ind,:);
dataft=abs(dataft).^2;
nuvec=nuvec(ind);


fitp=zeros(3,size(dataft,2));
fitp(1,:)=median(dataft,1);
fitp(2,:)=1.0;
fitp(3,:)=1.5;



for j=1:size(dataft,2);
  datasqr=dataft(:,j);
  fitp(:,j)=fminunc(@get_like_falpha,fitp(:,j));
end



function[like]=get_like_falpha(fitp,datasqr,nuvec)
amp=fitp(1);
knee=fitp(2);
slope=fitp(3);
global nuvec
global datasqr


model=amp*(1+(nuvec/knee).^(-slope));

like=sum(datasqr./model)+sum(log(model));

if ~isreal(like)
  like=1e10;
end
%like
