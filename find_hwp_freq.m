function[freqs]=find_hwp_freq(tods)
freqs=zeros(size(tods));
for j=1:length(tods),
  hwp=get_tod_hwp(tods);
  tvec=get_tod_tvec(tods);
  dt=diff(tvec);
  dhwp=diff(hwp);
  ii=find(dhwp>=0); %positive only, we'll just ignore wraps around zero here
  dhwp_tot=sum(dhwp(ii));
  dt_tot=sum(dt(ii));
  freqs(j)=dhwp_tot/dt_tot/(2*pi);
end
