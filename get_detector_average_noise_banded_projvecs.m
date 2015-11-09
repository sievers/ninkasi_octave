function[noise_per_samp]=get_detector_average_noise_banded_projvecs(tod)
%find the average white value for detectors from banded_projvecs noise model
%big_det_noise_amps(j,:)=mean(abs(myblock_clean).^2)/n;
%set_oneband_tod_noise_banded_projvec(tod,j,scale_facs(j)./big_det_noise_amps(j,:),vecs_use);
[whites,edges]=get_detector_white_noises_banded_projvec(tod);

%
wt_sum=zeros(size(whites,1),1);
for j=1:size(whites,2),
  wt_sum=wt_sum+whites(:,j)*(edges(j+1)-edges(j));
end
wt_sum=wt_sum/(edges(end)-1);
noise_per_samp=1./sqrt(wt_sum);

