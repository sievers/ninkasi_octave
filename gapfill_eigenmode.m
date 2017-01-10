function[dat]=gapfill_eigenmode(dat,tod,myopts)
thresh=get_struct_mem(myopts,'eiggap_thresh',4);
niter=get_struct_mem(myopts,'eiggap_niter',4);
nloop=get_struct_mem(myopts,'eiggap_nloop',3);

for loop=1:nloop,
  cc=dat'*dat;cc=0.5*(cc+cc');
  [vv,ee]=eig(cc);
  ee=diag(ee);
  ii=ee>(thresh*thresh*median(ee));
  vecs=vv(:,ii);
  mat=dat*vecs;
  %mat(:,end+1)=1;
  big_fitp=0;
  dat_tmp=dat;
  for j=1:niter, 
    fitp=inv(mat'*mat)*(mat'*dat_tmp);
    mdisp(['fit sum is ' num2str(sum(sum(abs(fitp))))])
    big_fitp=big_fitp+fitp;
    dat_tmp=dat_tmp-mat*fitp;
    
    push_tod_data(dat_tmp,tod);
    gapfill_data_c(tod);
    dat_tmp=get_tod_data(tod);
  end
  dat=dat_tmp+mat*big_fitp;
end


