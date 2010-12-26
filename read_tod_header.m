function[tod,lims,isrising]=read_tod_header(todname,pointing_name,decimate,varargin)
%interface to the C tod reader.

if ~exist('decimate')
  decimate=0;
end
if isempty(decimate)
  decimate=0;
end


assert(ischar(todname));
assert(~isempty(todname));
if ~isdir(todname)
  if ~exist([todname '.zip'])
    error(['cannot find tod ' todname ' or ' [todname '.zip'] ' in read_tod_header.'])
  end
end
%assert(isdir(todname));

if ~exist('pointing_name')
  pointing_name='';
end



point_tag=get_keyval_default('point_tag','2008',varargin{:});
point_dir=get_keyval_default('point_dir','/home/sievers/act/pointing/v06/',varargin{:});


if iscell('pointing_name')
  pointing_name=pointing_name{guess_tod_season(todname)};
end

if isstruct(pointing_name)
  tod=read_tod_header_nopoint_c(todname);
  myoffsets=get_tod_pointing_offsets_type_100(tod,pointing_name);
  set_tod_pointing_c(tod,myoffsets);
  
  return
end


if iscell(point_tag)
  point_tag=point_tag{guess_tod_season(todname)};
end

if iscell(point_dir)
  point_dir=point_dir{guess_tod_season(todname)};
end



if (isempty(pointing_name)),
  %pointing_name='/cita/d/raid-sievers/sievers/act/ninkasi/from_ishmael/pointing_offset_mbac145.txt';
  
  
  
  %pointing_name='/cita/d/raid-sievers/sievers/act/test_priors/obsfit_mbac145_south_rising_2008.for';
  [az,alt,ctime]=get_median_altaz(todname);
  [ra,dec]=alt_az_ctime2act_ra_dec_c(alt,az,ctime);
  [rah,ram,ras]=dd2dms(ra*180/pi/15);
  [dd,dm,ds]=dd2dms(dec*180/pi);
  


  if dec*180/pi>-25
    strip_equ='equ';
    isnorth=true;
  else
    strip_equ='south';
    isnorth=false;
  end

  
  isnorth=false;
  if (az>pi)
    if isnorth
      riseset='rising';
      isrising=true;
    else
      riseset='setting';
      isrising=false;      
    end
  else
    if isnorth
      riseset='setting';
      isrising=false;
      
    else
      riseset='rising';
      isrising=true;
    end
  end
  
  freq=get_tod_freq(todname);
  
  
  if (1)
    %pointing_name=['/home/sievers/act/pointing/v06/pointing_offset_mbac' freq '_' strip_equ  '_' riseset  '_2008.txt'];
    %pointing_name=['/home/sievers/act/pointing/v06/pointing_offset_mbac' freq '_' strip_equ  '_' riseset  '_' point_tag '.txt'];
    pointing_name=[point_dir '/pointing_offset_mbac' freq '_' strip_equ  '_' riseset  '_' point_tag '.txt'];
    mdisp(['apparent central position is ' rah ':' ram ':' ras '   ' dd ':' dm ':' ds '  using file ' pointing_name]);    
  else
    if (az>pi)
      
      %disp([todname ' is rising.']);
      
      %pointing_name='/home/sievers/act/pointing/pointing_offset_mbac145_south_rising_2008.txt';
      %pointing_name='/home/sievers/act/pointing/v06/pointing_offset_mbac145_south_rising_2008.txt';
      
      pointing_name=['/home/sievers/act/pointing/v06/pointing_offset_mbac' freq '_south_setting_2008.txt'];
      %pointing_name='/home/sievers/act/pointing/v06/pointing_nomdelt_mbac145_south_rising_2008.txt';
      
      
    else
      %disp([todname ' is setting.']);
      %pointing_name='/home/sievers/act/pointing/pointing_offset_mbac145_south_setting_2008.txt';
      %pointing_name='/home/sievers/act/pointing/v06/pointing_offset_mbac145_south_setting_2008.txt';
      
      %pointing_name='/home/sievers/act/pointing/v06/pointing_nomdelt_mbac145_south_setting_2008.txt';
      pointing_name=['/home/sievers/act/pointing/v06/pointing_offset_mbac' freq '_south_rising_2008.txt'];
      
      
      %pointing_name='/cita/d/raid-sievers/sievers/act/test_priors/v4/pointing_offset_mbac145_south_setting_2008.txt';
    end
  end
  
end

assert(ischar(pointing_name));



if (1)
  poff=read_pointing_offsets(pointing_name,varargin{:});
  if isempty(poff)
    poff=read_pointing_offsets('/home/sievers/act/pointing/v06_nomdelt_091215/pointing_offset_mbac145_south_setting_2008.txt');
    warning(['Unable to find offsets on file ' todname]);
    found_pointing=false;
  else
    found_pointing=true;
  end



  %mdisp('pointing offsets from structure.');
  if exist('decimate')
    [tod,lims]=read_tod_header_c(todname,poff,decimate);
  else
    [tod,lims]=read_tod_header_c(todname,poff);
  end


else
  if exist('decimate')
    [tod,lims]=read_tod_header_c(todname,pointing_name,decimate);
  else
    [tod,lims]=read_tod_header_c(todname,pointing_name);
  end
end

if ~found_pointing,
  lims=[];
end

createFFTWplans_c(tod);  %setup fft plans since we'll likely need
                         %them later.

