function[value]=cut_magic_carpets(tods,tol,min_cut)
%cut stationary ends of tods.  tol is in degrees, min_cut is in seconds.
if ~exist('tol')
  tol=0.2;
end

if ~exist('min_cut')
  min_cut=1;
end



tol=tol*pi/180;  %convert to radians

for j=1:length(tods),
  [alt,az]=get_tod_altaz(tods(j));

  nmin=ceil(min_cut/get_tod_dt(tods(j)));

  k=1;
  nsamp=length(az);
  while (abs(az(k)-az(1))<tol) & (k<nsamp)
    k=k+1;
  end
  if k>nmin,
    cut_tod_global_c(tods(j),1,k);
  end

  k=nsamp;
  while (abs(az(k)-az(nsamp))<tol) & (k>1)
    k=k-1;
  end

  if abs(k-nsamp)>nmin,
    cut_tod_global_c(tods(j),k,nsamp);
  end
end





  

  
