function[altaz]=get_detector_altaz(tod,row,col,exact)
if ~exist('col')
  det=row;
else
  [rr,cc]=get_tod_rowcol(tod);
  det=find((rr==row)&(cc==col));
end
disp(['working on detector ' num2str(det) '.']);
if isempty(det)
  error(['detector ' num2str(row) ',' num2str(col) ' does not exist in TOD.']);
else
  det=det-1;
end
if exist('exact')
  altaz=get_detector_altaz_c(tod,det,true);
else
  altaz=get_detector_altaz_c(tod,det);
end