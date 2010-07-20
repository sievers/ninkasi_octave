function[radec]=get_detector_radec(tod,row,col,exact)
if ~exist('col')
  det=row;
else
  [rr,cc]=get_tod_rowcol(tod);
  det=find((rr==row)&(cc==col));
end

if isempty(det)
  error(['detector ' num2str(row) ',' num2str(col) ' is not present in TOD.']);
else
  det=det-1;
end

if exist('exact')
  radec=get_detector_radec_c(tod,det,true);
else
  radec=get_detector_radec_c(tod,det);
end