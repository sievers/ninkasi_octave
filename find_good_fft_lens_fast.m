function[value]=find_good_fft_lens_fast(nmax,max_prime)
if ~exist('nmax','var')
  nmax=1001;
end
if ~exist('max_prime','var')
  max_prime=7;
end

vec=find_good_fft_lens_fast_c(nmax,max_prime);

value=sort(vec);

