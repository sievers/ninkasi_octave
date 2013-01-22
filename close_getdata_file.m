function[value]=close_getdata_file(df)
if ischar(df)
  return;
end

close_getdata_file_c(df);