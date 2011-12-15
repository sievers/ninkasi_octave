function[value]=apply_tod_noise_model(tod)

[myptr,mytype]=get_generic_tod_pointer(tod);

if mytype==0
  apply_tod_noise_model_c(tod);
  return
end
if (mytype==1)
  noise_mult_vis(myptr);
  return
end

