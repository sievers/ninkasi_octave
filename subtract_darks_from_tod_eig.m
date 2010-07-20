function[vecs,cc]=subtract_darks_from_tod_eig(tod,varargin)
%subtract off the dark detectors from the tod data.  Don't do anything fancy - just take all the dark 
%detector in a row and subtract off the average of 'em.

neig=get_keyval_default('neig_dark',7,varargin{:});
minfreq=get_keyval_default('minfreq_dark',5,varargin{:});


%JLS - 19Mar10, deal with potentially missing dark files
crud=get_dark_timestreams(tod);
if isempty(crud)
  vecs=[];
  cc=[];
  warning(['Dark subtraction failed due to missing file on ' get_tod_name(tod)]);
  return;
end

%darkdat=detrend_columns(get_dark_timestreams (tod));
darkdat=detrend_columns(crud);

%darkft_tmp=fft_r2c_octave(darkdat);
%darkft=fft_r2c_octave(darkdat);
darkft=fft_omp_r2c_many(darkdat);
%darkft=fft_r2c(darkdat);
do_even=iseven(size(darkdat,1));

total_t=size(darkdat,1)*get_tod_dt(tod);
nu=(1:size(darkft,1))'/total_t;
darkft(nu<minfreq,:)=0;
%dark=fft_c2r(darkft,do_even);clear darkft;
%dark=fft_c2r(darkft,do_even,true);clear darkft;
dark=fft_c2r_c(darkft,do_even);
amps=sqrt(mean(dark.^2,1));ind=(isfinite(amps)&(amps>0));amps=amps(ind);dark=dark(:,ind);dark=scale_rowcol_mat(dark,1./amps);

cc=dsyrk(dark);

if (neig>floor(0.5*length(cc)))
  neig_old=neig;
  neig=floor(0.5*length(cc));
  nm=get_tod_name(tod);
  if iscell(nm)
    nm=nm{1};
  end
  warning(sprintf('Trimming eigenvectors from %d to %d in darks on tod %s\n',neig_old,neig,nm));
end
[vv,ee]=eig(cc);
vv_use=vv(:,end-neig+1:end);vv_use=scale_rowcol_mat(vv_use,amps');

vecs=dark*vv_use;

data=detrend_columns(get_tod_data(tod));
if isempty(data)
  warning('tod is empty in subtract_darks_from_tod_blind');
  return;
end

%dataft=fft_r2c(data);clear data;
dataft=fft_r2c_octave(data);clear data;
dataft(nu<minfreq,:)=0;
%data=fft_c2r(dataft,do_even);clear dataft;
data=fft_c2r(dataft,do_even,true);clear dataft;

m1=vecs'*vecs;m2=vecs'*data;fitp=inv(m1)*m2;
clear data;

data=get_tod_data(tod)-vecs*fitp;
push_tod_data(data,tod);



