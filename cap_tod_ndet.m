function[ncut]=cap_tod_ndet(tod,ndet_max)
if ~exist('ndet_max')
  ndet_max=745;
end

[row,col]=get_tod_rowcol(tod);
ncut=0;
if length(row)>ndet_max,
  ncut=length(row)-ndet_max;
  for j=ndet_max+1:length(row),
    cut_detector_c(tod,row(j),col(j));
  end
  disp(['had too many detectors on ' get_tod_name(tod) ' cutting ' num2str(ncut) '.']);

  purge_cut_detectors(tod);
end