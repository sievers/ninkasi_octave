function[cal_facs,my_rr,my_ff]=get_tod_calib_factors(tod,varargin)
  fname=get_tod_name(tod); if fname(end)=='/', fname=fname(1:end-1);end;
%flatfield_dir=get_keyval_default('flatfield_dir','/cita/d/raid-sievers/sievers/act/shared/flatFields/2008v0/',varargin{:});
flatfield_dir=get_keyval_default('flatfield_dir','/home/sievers/act/shared/flatFields/2008v0/',varargin{:});
%response_dir=get_keyval_default('response_dir','/mnt/raid-cita/sievers/act/responsivities/2008/',varargin{:});
response_dir=get_keyval_default('response_dir','/home/sievers/act/responsivities/2008/',varargin{:});

mdisp(['getting flatfields from ' flatfield_dir]);
mdisp(['getting responsivities from ' response_dir]);

curdir=pwd;
if (fname(end-2:end)=='ar1'),
  cd(flatfield_dir);
  flatFieldAR1;
  cd(curdir);
else
  if (fname(end-2:end)=='ar2'),
    cd(flatfield_dir);
    flatFieldAR2;
    cd(curdir);
  else
    error('can''t find flatfields.');
  end
end



[myrows,mycols]=get_tod_rowcol(tod);
my_rr=zeros(length(myrows),1);
my_ff=zeros(length(myrows),1);


[a,b,c]=fileparts(fname);
response_name=[response_dir '/'  b  c];
try

  responsivity=load(response_name);

  %responsivity=load([response_dir '/'  b  c]);


  for j=1:length(myrows),
    my_rr(j)=responsivity( (responsivity(:,1)==myrows(j))&(responsivity(:,2)==mycols(j)),3);
    ind=(col==mycols(j))&(row==myrows(j));
    if sum(ind)==0
      my_ff(j)=nan;
    else
      my_ff(j)=cal(ind);
    end
    
  end
catch
  my_rr(:)=nan;
  my_ff(:)=nan;
  warning(['responsivity file ' response_name ' not found']);
end

mdisp(['Flatfield: ' flatfield_dir '  Responsivity: ' response_dir]);

cal_facs=my_rr.*my_ff;
if sum(isnan(cal_facs))>0,
  mdisp(['have uncalibrated detectors that aren''t cut.  Cutting now on ' fname '.']);
  for j=1:length(cal_facs),
    if isnan(cal_facs(j)),
      cut_detector_c(tod,myrows(j),mycols(j));
    end
  end
end



  
  
