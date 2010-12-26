function[data_filt,filtered,big_minv]=get_timestream_chisq_modeproj(data,data_clean,vecs,amps,to_filt)
%function[data_filt,filtered,big_minv]=get_timestream_chisq_modeproj(data,data_clean,vecs,amps,to_filt)
%vecs is, for instance, our estimate of the common mode.  Note that this is *only* used for 
%noise estimation, not for actual removal from timestreams.  Amps is the response of each detector to vecs.
%to_filt is a cell array containing the model timestreams one wishes to apply ninv to.


%make common mode truly common
%amps=ones(size(amps))*mean(mean(amps));

dd=abs(fft(data_clean).^2);
nn=floor((size(dd,1))/2+1);
dd=dd(1:nn,:);
nmin=round(0.05*size(dd,1));
if nmin<2,
  nmin=2;
end


sigsqr=mean(dd(2:end,:),1);
%make variances equal:
%sigsqr=ones(size(sigsqr))*mean(sigsqr);

clear dd
clear data_clean;
wt=1./sigsqr;

tic
smooth_ps=get_logsmoothed_ps(vecs);

datft=fft(data);
datft=datft(1:nn,:);
datft_filt=datft;
datft_filt(1,:)=0;
ampsqr=amps.^2;
a=diag(wt);
ainv=diag(1./wt);
au=amps*ainv;
noise=diag(sigsqr);

for j=1:length(to_filt),
  crap=fft(to_filt{j});
  crapft=crap(1:nn,:);
  crapft(1,:)=0;
  to_filt_ft(j)={crapft};
end

%let's just brute-force it here, should be quick enough for tens of detectors
for j=length(smooth_ps):-1:2,
  mat=amps'*diag(smooth_ps(j,:))*amps;
  minv=inv(mat+noise);
  datft_filt(j,:)=datft_filt(j,:)*minv;
  for k=1:length(to_filt_ft),
    to_filt_ft{k}(j,:)=to_filt_ft{k}(j,:)*minv;
  end
end
toc


if iseven(size(data,1)),
  crud=datft_filt;crud=[crud;flipud(conj(datft_filt(2:end-1,:)))];data_filt=real(ifft(crud));
  for j=1:length(to_filt_ft),
    crud=to_filt_ft{j};crud=[crud;flipud(conj(crud(2:end-1,:)))];filtered(j)={real(ifft(crud))};
  end
else
  crud=datft_filt;crud=[crud;flipud(conj(datft_filt(2:end,:)))];data_filt=real(ifft(crud));
  for j=1:length(to_filt_ft),
    crud=to_filt_ft{j};crud=[crud;flipud(conj(crud(2:end,:)))];filtered(j)={real(ifft(crud))};
  end
end
