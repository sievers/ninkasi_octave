function[meds]=omp_median(dat,dim)
%broken right now, use this to get around it.
if (1)
  meds=median(dat,dim);
  return
end
disp('hi')

if ndims(dat)>2
  error('cannot deal with more than 2 dimensions in omp_median.');
end


if ~exist('dim'),
  if size(dat,1)==1,
    meds=omp_median_r(dat);
    return;
  end
  if size(dat,2)==1,
    meds=omp_median_c(dat);
    return;
  end
  dim=1;
end

if (dim<1) |(dim>2)
  error('dim must be 1 or 2');
end
if (dim==1)
  meds=omp_median_c(dat);
else
  meds=omp_median_r(dat);
end


