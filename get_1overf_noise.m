function[value,spec,nu]=get_1overf_noise(ndata,ndet,white,amp_1,ind,dt,pad)
if ~exist('pad')
    pad=false;
end

if ~exist('dt')
    dt=1/400;
end
if ~exist('ind')
    ind=-1;
end
if ~exist('white')
    white=1;
end
if ~exist('amp_1')
    amp_1=white;
end


if pad,
    ndata=ndata*2;
end


if iseven(ndata),
    nu=fftshift(-ndata/2:ndata/2-1)';
else
    nu=ifftshift( ((-ndata+1)/2):((ndata-1)/2))';
end
nu=abs(nu)/(ndata*dt);
nu(1)=0.5*nu(2);  %get a non-zero value in here

tic
crud=fft2(randn(ndata,ndet));
toc

spec=white^2+  amp_1^2*((nu/1).^ind);
if exist('lowmem')
    for j=1:ndet,
        crud(:,j)=crud(:,j).*sqrt(spec);
    end
else
    crud=crud.*repmat(sqrt(spec),[1 ndet]);
end
tic
value=real(ifft2(crud));
toc
if pad,
    value=value(1:end/2,:);
    spec=spec(1:2:end);
    nu=nu(1:2:end);
end




