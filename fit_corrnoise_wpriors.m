function[fit,mapft,fitft]=fit_corrnoise_wpriors(vecs,map,prior)
little_mat=vecs*vecs';
mapft=map;
for j=1:size(map,2),
  mapft(:,j)=fft(mapft(:,j));
end
if (~exist('prior'))
  fitft=mapft*inv(little_mat);
else
  val1=multiply_corrnoise_wprior_c(real(mapft)',prior',little_mat);
  val2=multiply_corrnoise_wprior_c(imag(mapft)',prior',little_mat);
  fitft=val1'+i*val2';
  %for j=1:length(map),
  %  fitft(j,:)=mapft(j,:)*inv(little_mat+diag(1./prior(j,:).^2));
  %end
end

fit=zeros(size(fitft));
for j=1:size(fitft,2),
  fit(:,j)=real(ifft(fitft(:,j)));
end

  

