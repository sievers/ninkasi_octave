function[value]=fft_c2r(mat,even,do_new)
if ~exist('even')
  even=isreal(mat(end,:));
end
if ~exist('do_new')
  do_new=false;
end

if (0)  %do transform locally

  n=size(mat,1);
  m=round(size(mat,2)/2);
  %do things in half-blocks so that memory doesn't overflow
  value1=mat(:,1:m);
  value2=mat(:,m+1:end);
  if even,
    value1(end+1:end+n-2,:)=flipud(conj(mat(2:end-1,1:m)));
    value2(end+1:end+n-2,:)=flipud(conj(mat(2:end-1,m+1:end)));
  else
    value1(n+1:n+n-1,:)=flipud(conj(mat(2:end,1:m)));
    value2(n+1:n+n-1,:)=flipud(conj(mat(2:end,m+1:end)));
  end
  value=[real(ifft(value1)) real(ifft(value2))];
else

  if even
    nn=2*size(mat,1)-2;
  else
    nn=2*size(mat,1)-1;
  end
  value=fft_c2r_c(mat,even)/nn;
end  
