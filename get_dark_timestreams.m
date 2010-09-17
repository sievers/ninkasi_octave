function[dat,rows,cols]=get_dark_timestreams(tod,dirroot)
if exist('dirroot')
  [rows,cols]=get_darkdets_rowcol(tod,dirroot);
else
  [rows,cols]=get_darkdets_rowcol(tod);
end


if (~isempty(rows))
  dat=read_dirfile_tod_data_from_rowcol_list_c(tod,rows,cols);
else
  dat=[];
end


