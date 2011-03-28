function[ncut]=cap_tod_ndet(tod,ndet_max)
if numel(tod)>1,
  for j=1:length(tod),
    if exist('ndet_max')
      cap_tod_ndet(tod(j),ndet_max);
    else
      cap_tod_ndet(tod(j));
    end
  end
  return
end


if ~exist('ndet_max')
  if (1)
    %matrices don't want to get larger than 2^31 8-byte elements
    ndata=get_tod_ndata(tod);
    ndet_max=floor(2^28/(ndata+2));
  else
    ndet_max=745;
  end
end

[row,col]=get_tod_rowcol(tod);
ncut=0;
if length(row)>ndet_max,
  ncut=length(row)-ndet_max;
  for j=ndet_max+1:length(row),
    cut_detector_c(tod,row(j),col(j));
  end
  disp(['had too many detectors on ' get_tod_name(tod) ' cutting ' num2str(ncut) ', trimming to ' num2str(ndet_max) '.']);

  purge_cut_detectors(tod);
end