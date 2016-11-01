function[sin_fitp,cos_fitp]=set_tod_twogamma_fit(tod,varargin)
if (does_tod_have_twogamma_fit(tod))
  mdisp('tod already has twogamma_fit in set_tod_twogamma_fit');
  return
end
downsamp=get_keyval_default('downsamp',0,varargin{:}); %how tightly we should subsample
npoly=get_keyval_default('npoly_2gamma',3,varargin{:});
rescale_az=get_keyval_default('rescale_az_2gamma',false,varargin{:});

if (downsamp>1),
  aaa=now;
  [myra,mydec,twogamma]=precalc_actpol_pointing_exact_subsampled_c(tod,downsamp);
  bbb=now;
  clear myra;
  clear mydec;
  [alt,az]=get_tod_altaz(tod);
  ndata=get_tod_ndata(tod);
  tvec=(0:(ndata-1))'/ndata;

  %tvec=tvec(1:downsamp:(end-downsamp+1));
  %az=az(1:downsamp:(end-downsamp+1));
  %alt=alt(1:downsamp:(end-downsamp+1));


  tvec=tvec(1:downsamp:end);
  az=az(1:downsamp:end);
  if (rescale_az)
    az_cent=mean(az);
    az=az-az_cent;
    az_std=std(az);
    az=az/az_std;
  else
    az_cent=0;
    az_std=1;
  end


  alt=alt(1:downsamp:end);

  mat=ones(length(tvec),2+npoly);
  mat(:,1)=1;
  for j=1:npoly,
    mat(:,j+1)=mat(:,j).*az;
  end
  mat(:,end)=tvec;
  mm=mat'*mat;
  sin_rhs=mat'*sin(twogamma);
  sin_fitp=invscale(mm)*sin_rhs;
  clear sin_rhs;
  cos_rhs=mat'*cos(twogamma);
  cos_fitp=invscale(mm)*cos_rhs;
  clear cos_rhs
  ccc=now;
  disp(86400*([bbb-aaa ccc-bbb]));
  set_tod_twogamma_fit_c(tod,sin_fitp,cos_fitp,az_cent,az_std);
  return
end





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

