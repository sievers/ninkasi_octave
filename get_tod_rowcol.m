function[rows,cols]=get_tod_rowcol(tod)
assert(class(tod)=='int64');
rows=get_tod_row_c(tod);
cols=get_tod_col_c(tod);
