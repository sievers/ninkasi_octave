function[asdfadf]=set_tod_twogamma_fit(tod,varargin)
if (does_tod_have_twogamma_fit(tod))
  disp('tod already has twogamma_fit in set_tod_twogamma_fit');
  return
end


npoly=get_keyval_default('npoly_2gamma',3,varargin{:});

try
  twogamma=get_tod_2gamma(tod);
catch
  mdisp('calculating pointing inside set_tod_twogamma_fit.')
  precalc_two_gamma=2;  %flag set in ninkasi
  precalc_actpol_pointing_exact(tod,precalc_two_gamma);
  twogamma=get_tod_2gamma(tod);
  free_tod_2gamma_saved(tod);  
end
[alt,az]=get_tod_altaz(tod);
ndata=get_tod_ndata(tod);
tvec=(0:(ndata-1))'/ndata;
mat=ones(ndata,2+npoly);
mat(:,1)=1;
for j=1:npoly,
  mat(:,j+1)=mat(:,j).*az;
end
mat(:,end)=tvec;
mm=mat'*mat;
sin_rhs=mat'*sin(twogamma);
sin_fitp=inv(mm)*sin_rhs;
clear sin_rhs;
cos_rhs=mat'*cos(twogamma);
cos_fitp=inv(mm)*cos_rhs;
clear cos_rhs
if (0)
  max_sinerr=max(max(mat*sin_fitp-sin(twogamma)));
  max_coserr=max(max(mat*cos_fitp-cos(twogamma)));
  mdisp(['Max cos/sin errors in set_tod_twogamma_fit are ' num2str([max_coserr max_sinerr])]);
  clear twogamma
end

%disp([sin_fitp(:,1) cos_fitp(:,1)])
set_tod_twogamma_fit_c(tod,sin_fitp,cos_fitp);

