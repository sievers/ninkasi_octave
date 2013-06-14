function[val]=invert_1_over_f_common_mode_data(data,knees,whites,inds,cm_vec,cm_white,cm_knee)
tic;
n=size(data,1);
if iseven(n)
  nn=(n/2)+1;
else
  nn=(n+1)/2;
end
cm_vec=rowvec(cm_vec);
toc
datft=fft(data);
toc
datft=datft(1:nn,:);
nu=(1:nn)'-1;nu(1)=0.1;  %make sure it's not entirely singular
if size(whites,2)==1,
  whites=whites';
end
amat=repmat(whites.^-2,[nn 1]);
for j=1:size(amat,2),
  amat(:,j)=amat(:,j)./(1+  (nu/knees(j)).^inds(j));
end
amat(1,:)=0;  %explicitly chop the mean of the timestream
toc
ninv_datft=zeros(size(datft));

%clear;n=5;a=diag(randn(n,1));x=randn(n,1);a=a+a';v=randn(n,1);ainv=inv(a);ai=diag(ainv);c=sum(v.^2.*ai);rhs=[x x+5*v x+1e3*v];y=rhs.*repmat(ai,[1 3]);yy=sum(repmat(v,[1 3]).*y)/c;yy2=y-ainv*(v*yy)
%c=sum(v.^2.*ai);rhs=x;y=rhs.*ai;yy=sum(v.*y)/c;yy2=y-ai.*v.*yy
if ~exist('cm_white')
  for j=2:nn,
    ai=amat(j,:);
    %make a^-1 x first
    y=datft(j,:).*ai;
    c=sum(cm_vec.^2.*ai);
    yy=sum(cm_vec.*y)/c;
    ninv_datft(j,:)=y-ai.*cm_vec.*yy;
  end
end
toc

if iseven(n)
  ninv_datft=[ninv_datft;flipud(conj(ninv_datft(2:end-1,:)))];
else
  ninv_datft=[ninv_datft;flipud(conj(ninv_datft(2:end,:)))];
end
val=real(ifft(ninv_datft));
toc 
