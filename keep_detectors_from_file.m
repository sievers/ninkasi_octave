function[value]=keep_detectors_from_file(tods,fname)
lines=read_text_file_comments(fname);
keeprows=zeros(numel(lines),1);
keepcols=keeprows;
for j=1:numel(lines),
  crud=sscanf(lines{j},'%d');
  assert(numel(crud)==2);
  keeprows(j)=crud(1);
  keepcols(j)=crud(2);
end
keepvec=keeprows+i*keepcols;
for j=1:numel(tods),
  [rows,cols]=get_tod_rowcol(tods(j));
  for k=1:numel(rows),
    vv=rows(k)+i*cols(k);
    if min(abs(keepvec-vv))>0
      cut_detector_c(tods(j),rows(k),cols(k));
    end
  end
end
