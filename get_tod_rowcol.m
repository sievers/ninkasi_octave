function[rows,cols]=get_tod_rowcol(tod)
try
  assert(class(tod)=='int64');
catch
  disp(['bad class of tod: ' class(tod)])
  assert(1==0);
end

rows=get_tod_row_c(tod);
cols=get_tod_col_c(tod);
