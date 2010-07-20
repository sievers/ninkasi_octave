function[mat]=fft_r2c_octave(mat)
if (1)
  mat=fft_omp_r2c_many(mat);
else
  m=round(size(mat,2)/2);

  if iseven(size(mat,1))
    nn=1+size(mat,1)/2;
  else
    nn=(1+size(mat,1))/2;
  end
  mat1=fft(mat(:,1:m));
  mat1=mat1(1:nn,:);
  mat=fft(mat(:,m+1:end));
  mat=mat(1:nn,:);
  mat=[mat1 mat];
end
