function[asdf]=fit_hwp_az_poly_to_data(tod,myopts)
niter=get_struct_mem(myopts,'hwp_niter',3);
npoly=get_struct_mem(myopts,'hwp_npoly',4);
nsin=get_struct_mem(myopts,'hwp_nsin',30);
naz=get_struct_mem(myopts,'hwp_naz',2);

for j=1:niter,
  gapfill_data_c(tod);
  fit_hwp_az_poly_to_data_c(tod,nsin,naz,npoly);
end
gapfill_data_c(tod);




