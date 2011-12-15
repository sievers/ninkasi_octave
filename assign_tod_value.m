function[value]=assign_tod_value(tod,val)
[myptr,mytype]=get_generic_tod_pointer(tod);
if mytype==0
  assign_tod_value_c(tod,val);
  return
end
if (mytype==1)
  assign_vis_value(myptr,val);
  return
end


