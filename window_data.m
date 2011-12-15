function[value]=window_data(tods)
for j=1:length(tods),
  mytype=query_tod_type(tods(j));
  if (mytype==0)
    window_data_c(tods(j));
  end
end
