function[value]=allocate_tod_storage(tod)
mytype=query_tod_type(tod);
if mytype==0
  allocate_tod_storage_c(tod);
end

