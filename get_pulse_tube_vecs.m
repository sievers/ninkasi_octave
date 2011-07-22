function[vecs]=get_pulse_tube_vecs(tod,freq,nufrac,nvec)
if ~exist('freq')
  freq=[];
end
if isempty(freq)
  freq=1.411;
end
if ~exist('nufrac')
  nufrac=[];
end
if isempty(nufrac)
  nufrac=0.4;
end
if ~exist('nvec')
  nvec=[];
end
if isempty(nvec)
  nvec=1;
end


tvec=get_tod_tvec(tod);


nuvec=get_tod_nuvec(tod);dnu=nufrac*(nuvec(3)-nuvec(2));
%dnu=1/(tvec(end)-tvec(1))*nufrac;
freqs=freq+(-nvec:nvec)*dnu
nn=numel(freqs);
for j=nn:-1:1,
  phi(:,j)=2*pi*tvec*freqs(j);
end
vecs=[cos(phi) sin(phi)];