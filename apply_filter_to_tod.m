function[value]=apply_filter_to_tod(tod,filt)
n=get_tod_ndata(tod);
assert(length(filt)> (n/2));

if isreal(filt)
  apply_real_filter_to_data(tod,filt);
else
  apply_complex_filter_to_data(tod,filt);
end
