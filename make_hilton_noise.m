function[newmat]=make_hilton_noise(mat,varargin)

newmat=0;
nmode=get_keyval_default('nmode',12,varargin{:});
mdisp(['want ' num2str(nmode) ' modes in make_hilton_noise.'])
if nmode>size(mat,2)
  warning('too many modes requested in make_hilton_noise.  Truncating.');
  nmode=size(mat,2);
end


mywhites=omp_median(abs(fft_r2c(mat)),1)/sqrt(size(mat,1)/sqrt(2));

if (0)
  if sum(sum(~isfinite(mywhites)))
    warning('mywhites has nastiness in make_hilton_noise.')
    newmat=0;
    return;
  end
end

if nmode>0
  %mycorr=mat'*mat;  mycorr=mycorr+mycorr';


  if sum(sum(~isfinite(mat)))
    warning('mat has nastiness in make_hilton_noise')
    newmat=0;
    return;
  end


  mycorr=dsyrk(mat);
  if sum(sum(~isfinite(mycorr)))
    warning('mycorr has nastiness in make_hilton_noise')
    newmat=0;
    whos
    fid=fopen('nasty.dat','w');
    fwrite(fid,size(mat),'int');
    fwrite(fid,mat,'double');
    fclose(fid);
    return;
  end
  [v,e]=eig(mycorr);
  vv=v(:,end-nmode+1:end);
  mymodes=mat*vv;
  mymodes_ft=fft_r2c(mymodes);
  signvec=randn(size(mymodes_ft,1),1);signvec=2*(signvec>0)-1;
  mymodes_ft=scale_rowcol_mat(real(mymodes_ft),signvec)+i*scale_rowcol_mat(imag(mymodes_ft),signvec);
  mymodes=fft_c2r(mymodes_ft,iseven(size(mat,1)));  %OK, got the phase-scattered modes.
  newmat=scale_rowcol_mat(randn(size(mat)),mywhites)+mymodes*vv';
else
  newmat=scale_rowcol_mat(randn(size(mat)),mywhites);
end
