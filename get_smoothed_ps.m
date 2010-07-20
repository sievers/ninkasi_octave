function[ps]=get_smoothed_ps(vec,len)
if ~exist('len')
  len=20;
end

if (size(vec,2)>1),
  ps=0*vec;
  for j=1:size(vec,2),
    ps(:,j)=get_smoothed_ps(vec(:,j),len);
  end
  return
end

    
vecft=abs(fft(vec));
vecft=[vecft;vecft;vecft];
vecft=10.^smooth_data_fft(log10(vecft),len);
ps=real(vecft(end/3+1:2*end/3));

