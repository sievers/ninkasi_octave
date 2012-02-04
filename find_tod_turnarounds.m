function[turns]=find_tod_turnarounds(tod)
[alt,az]=get_tod_altaz(tod);
az_smooth=real(smooth_data_fft(az,30));

crap1=diff(az_smooth(1:end-1));crap2=diff(az_smooth(2:end));crud=crap1.*crap2;

ind=1:length(crud);fwee=crud<0;
tol=0.01;
fwee=fwee&((az(2:end-1)>max(az)-tol)|(az(2:end-1)<min(az)+tol));

turns=ind(fwee);
flub=unique(diff(turns));
if (std(flub)>3)
  warning(['Probably had some failures in find_tod_turnarounds.  Scan lengths in samples were ' num2str(flub) ' in file ' get_tod_name(tod) ]);
end


