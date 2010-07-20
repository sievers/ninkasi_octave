function[value]=select_detector_columns(tod,columns)
[myrows,mycols]=get_tod_rowcol(tod);
ncut=0;
for j=1:length(mycols),
  if sum(mycols(j)==columns)==0
    ncut=ncut+1;
    cut_detector_c(tod,myrows(j),mycols(j));
  end
end
mdisp(['cut ' num2str(ncut) ' detectors in select_detector_columns.']);
