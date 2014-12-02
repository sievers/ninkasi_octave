function[tod,isok]=read_tod_header_actpol(tod_name,varargin)
split_cells=get_keyval_default('split_cell_names',true,varargin{:});
if (iscell(tod_name)&(split_cells)),
  for j=1:length(tod_name),
    [tod(j),isok(j)]=read_tod_header_actpol(tod_name{j},varargin{:});
  end
  return
end
tod_name=_strip_zips(tod_name);

%/project/r/rbond/sievers/actpol/calib/actpol1_2013_c1_cal1_atm
calib_dir=get_keyval_default('calib_dir','',varargin{:});
calib_tail=get_keyval_default('calib_tail','calib.txt',varargin{:});
tod_offsets=get_keyval_default('tod_offsets',[],varargin{:});
offsets=get_keyval_default('offsets',[],varargin{:});
cuts=get_keyval_default('cuts','',varargin{:});
if isempty(offsets)
  error(['must have offsets defined']);
end
if size(offsets,2)<4
  error(['must have row,col,dx,dy in offsets']);
end


%try

  %if we don't have polarization info, set angles to zero.
  if size(offsets,2)==4,
    offsets(:,end+1)=0;
  end


  rr=offsets(:,1);
  cc=offsets(:,2);
  dx=offsets(:,3);
  dy=offsets(:,4);
  theta=offsets(:,5);

  [cpu_s,cpu_us,az_raw,el_raw,flags]=read_many_dirfile_channels(tod_name,'cpu_s','cpu_us','az','el','enc_flags');
  if ~isempty(cuts)
    [nsamp,samp_offset]=read_cuts_octave([],cuts,varargin{:})

    
    last_samp=samp_offset+nsamp;
    if last_samp>numel(cpu_s)
      warning('too many samples coming out of cuts.  Truncating.')
      last_samp=numel(cpu_s);
    end
    %this block should be good to use, leave bottom just in case
    if (1)
      cpu_s=cpu_s(samp_offset+1:last_samp);
      cpu_us=cpu_us(samp_offset+1:last_samp);
      az_raw=az_raw(samp_offset+1:last_samp);
      el_raw=el_raw(samp_offset+1:last_samp);
      flags=flags(samp_offset+1:last_samp);     
    else
      cpu_s=cpu_s(samp_offset+1:samp_offset+nsamp);
      cpu_us=cpu_us(samp_offset+1:samp_offset+nsamp);
      az_raw=az_raw(samp_offset+1:samp_offset+nsamp);
      el_raw=el_raw(samp_offset+1:samp_offset+nsamp);
      flags=flags(samp_offset+1:samp_offset+nsamp);
    end
  else
    samp_offset=0;
  end

  ct=cpu_s+1e-6*cpu_us;
  if sum(ct==0)>0
    ii=find(ct==0);
    ct=_repair_ct(ct);
    %the samples with ctime=0 seem to also have az=0 and el=0
    %and they don't seem to be getting picked up by the flags
    az_raw(ii)=median(az_raw);
    el_raw(ii)=median(el_raw);
  end

  az=az_raw*2.682209174765e-06 - 1.825699070000e+01;
  el=el_raw*-2.682209174765e-06 + 2.119976743000e+02;
  badinds=find(((flags==8)|(flags==24)));

  for jj=1:length(badinds),
    if badinds(jj)==1,
      az(badinds(jj))=median(az);
      el(badinds(jj))=median(el);
    else
      az(badinds(jj))=az(badinds(jj)-1);
      el(badinds(jj))=el(badinds(jj)-1);
    end
  end



  tod=allocate_tod_c();
  
  
  set_tod_ndata_c(tod,length(az));
  set_tod_nsamp_offset_c(tod,samp_offset);
  set_tod_altaz_c(tod,el*pi/180,az*pi/180);
  set_tod_timevec_c(tod,ct);
  set_tod_dt_c(tod,median(diff(ct)));
  set_tod_rowcol_c(tod,rr,cc);
  set_tod_filename(tod,tod_name);
  if iscell(tod_name), tod_name=tod_name{1};end;
  alloc_tod_cuts_c(tod);
  
  %tod_offsets=read_text_file_comments(['/home/mhasse/depot/TODOffsets/' field_tag '_130916/offsets.txt']);
  %tod_offsets=parse_tod_offsets(tod_offsets);
  %[tod_dx,tod_dy]=get_tod_offset(tod_names,tod_offsets);


  if ~isempty(tod_offsets)
    [tod_dx,tod_dy]=get_tod_offset(tod_name,tod_offsets);
    if isnan(tod_dx)|isnan(tod_dy)
      warning(['offsets not found on ' tod_name]);
      isok=false;
      return
    end
  else
    tod_dx=0;
    tod_dy=0;
  end
  initialize_actpol_pointing(tod,dy+tod_dy,dx+tod_dx,theta,148.0,1);

  if ~isempty(calib_dir)
    tod_name
    tt=tod_name;
    ii=find(tt=='/');
    if ~isempty(ii)
      if ii(end)==length(tt),
        tt=tt(1:end-1);
        ii=ii(1:end-1);
      end
      tt=tt(ii(end)+1:end);
    end
    %calfile=[calib_dir '/' tt '/calib.txt'];
    %calfile=[calib_dir '/' tt '.calib.txt'];
    calfile=[calib_dir '/' tt '.' calib_tail];

    cals=read_python_dict(calfile);
    cals.rr=floor(cals.det_uid/32);
    cals.cc=cals.det_uid-32*cals.rr;
    nr=max([max(cals.rr) max(rr)])+1;
    nc=max([max(cals.cc) max(cc)])+1;
    cal_mat=zeros(nr,nc);
    for j=1:length(cals.rr),
      cal_mat(cals.rr(j)+1,cals.cc(j)+1)=cals.cal(j);
    end
    det_resp=0*rr;
    for j=1:length(rr),
      det_resp(j)=cal_mat(rr(j)+1,cc(j)+1);
    end
    save_calib_facs_c(tod,det_resp);
    purge_cut_detectors(tod);
    
  end
  isok=true;
%catch
if (0)
  warning(['tod ' tod_name ' failed in read_tod_header_actpol for some reason.']);
  if exist('tod')
    destroy_tod(tod);
  end
  isok=false;
end



function[tod_name]=_strip_zips(tod_name)
tag='.zip';
nn=numel(tag);
if ischar(tod_name)
  if strcmp(tod_name(end-nn+1:end),tag)
    tod_name=tod_name(1:end-nn);
    disp(tod_name);
  end
end

function[vec]=_repair_ct(vec)
ii=find(vec==0);
dt=median(diff(vec));
assert(ii(1)>1);
for j=1:length(ii),
  vec(ii(j))=vec(ii(j)-1)+dt;
end

