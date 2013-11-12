function[value]=close_getdata_file(df)
if ischar(df)
  return;
end

for j=1:length(df),
  close_getdata_file_c(df(j));
end
