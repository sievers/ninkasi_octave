function[value]=cut_column(tods,cols_in)

if length(tods)>1,
  for j=1:length(tods),
    cut_column(tods(j),cols_in);
  end
  return
end


[rows,cols]=get_tod_rowcol(tods);

for k=1:length(cols_in),
  for j=1:length(rows),
    if cols(j)==cols_in(k)
      cut_detector_c(tods,rows(j),cols(j));
    end
  end
end
