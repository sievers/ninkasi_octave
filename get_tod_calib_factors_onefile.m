function[cal]=get_tod_calib_factors_onefile(tod,mydir,opts)

tail_tag=get_struct_mem(opts,'calib_tail','cal0');
tod_name=get_tod_name(tod);
tags=strsplit(strtrim(tod_name),'/',true);
tag=tags{end};
if isempty(tag)
  tag=tags{end-1};
end

fname=[mydir '/' tag '.' tail_tag];
mdisp(['getting new calibration factors from ' fname]);
[cal_in,row_in,col_in]=read_calib_factors(fname);
if isempty(cal_in)
  %error(['unable to find file ' fname ' in get_tod_calib_factors_onefile.']);
end
[rows,cols]=get_tod_rowcol(tod);
cal=0*rows;cal(1:numel(cal))=nan;

%assign calib factors to found detectors
for j=1:length(row_in),
  ind=find((rows==row_in(j))&(cols==col_in(j)));
  if ~isempty(ind),
    cal(ind)=cal_in(j);
  end
end

if sum(isnan(cal)>0)
  mdisp(['cutting ' num2str(sum(isnan(cal))) ' uncalibrated detectors on ' tag]);
  for j=1:length(cal)
    if isnan(cal(j))
      cut_detector_c(tod,rows(j),cols(j));
    end
  end
end
