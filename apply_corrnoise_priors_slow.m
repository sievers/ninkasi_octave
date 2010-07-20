function[fit,mapft,fitft]=apply_corrnoise_priors_slow(vecs,map,prior)
little_mat=vecs*vecs';
mapft=map;
for j=1:size(map,2),
  mapft(:,j)=fft(mapft(:,j));
end
if (~exist('prior'))
  fitft=mapft*inv(little_mat);
else
  fitft=0*mapft;
  for j=1:length(mapft),
    fitft(j,:)=mapft(j,:)*inv(little_mat+diag(1./prior(j,:).^2));
  end  
end

fit=zeros(size(fitft));
for j=1:size(fitft,2),
  fit(:,j)=real(ifft(fitft(:,j)));
end

  

