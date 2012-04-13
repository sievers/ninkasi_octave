function[value]=get_alldata_fft_plans(tods,wisdomname)
if ~exist('wisdomname')
  %wisdomname='/scratch/sievers/.fftw_wisdom';
  wisdomname='/project/r/rbond/sievers/.fftw_wisdom';
end
for j=1:length(tods),
  if ~isempty(wisdomname)
    get_alldata_fft_plans_c(tods(j),wisdomname)
  else
    get_alldata_fft_plans_c(tods(j));
  end
end
