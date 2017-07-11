function [dets, tods]=read_tod_spider(dets, samp_min, samp_max, varargin)

% Default options
data_root = get_keyval_default('data_root', getenv('SPIDER_DATA_ROOT'), varargin{:});
spider_tools_path = get_keyval_default('spider_tools_path', getenv('SPIDER_TOOLS_PATH'), varargin{:});
ver_point = get_keyval_default('ver_point', 'point03', varargin{:});
ver_time = get_keyval_default('ver_time', 'time02', varargin{:});
ver_clean = get_keyval_default('ver_clean', 'dcclean08', varargin{:});
ver_flag_comb = get_keyval_default('ver_flag_comb', 'flag_comb08', varargin{:});
ver_flag_det_fj = get_keyval_default('ver_flag_det_fj', 'flag_det_fj06', varargin{:});
ver_flag_det_jump = get_keyval_default('ver_flag_det_jump', 'flag_det_jump06', varargin{:});
ver_gain = get_keyval_default('ver_gain', 'gain09', varargin{:});
ver_hwp = get_keyval_default('ver_hwp', 'hwp03', varargin{:});
gain_tag = get_keyval_default('gain_tag',ver_gain,varargin{:});
if isempty(gain_tag)
  gain_tag=ver_gain;
end

centroid_file = get_keyval_default('centroid_file', [spider_tools_path '/analysis/config/bolotables/centroids.txt'], varargin{:});
angle_file = get_keyval_default('angle_file', [spider_tools_path '/analysis/config/bolotables/trpns_pol_angles_measured.txt'], varargin{:});
jump_bits = get_keyval_default('jump_bits', [16 17], varargin{:});
min_frac = get_keyval_default('min_frac', 0.5, varargin{:});
poly_order = get_keyval_default('poly_order', 20, varargin{:});
tod_filename = get_keyval_default('tod_filename', 'spider', varargin{:});
%-------------------------------------------------------

ndet=size(dets,1);

% Initialize dirfiles
aaa=now;
pointfile = init_getdata_file([data_root '/' ver_point]);
datfile = init_getdata_file([data_root '/' ver_clean]);
timefile = init_getdata_file([data_root '/' ver_time]);
flagfile = init_getdata_file([data_root '/' ver_flag_comb]);
fjumpfile = init_getdata_file([data_root '/' ver_flag_det_fj]);
jumpfile = init_getdata_file([data_root '/' ver_flag_det_jump]);
gainfile = init_getdata_file([data_root '/' ver_gain]);
hwpfile = init_getdata_file([data_root '/' ver_hwp]);
bbb=now;
disp(['took ' num2str((bbb-aaa)*86400) ' seconds to initialize files'])


% Read pointing
ra_bore = getdata_double_channel_piece(pointfile,['RA_' upper(ver_point)],'first_sample',samp_min,'num_samples',samp_max-samp_min);
ra_bore=ra_bore*15;  % Convert RA from hours to degrees
dec_bore = getdata_double_channel_piece(pointfile,['DEC_' upper(ver_point)],'first_sample',samp_min,'num_samples',samp_max-samp_min);
phi_bore = getdata_double_channel_piece(pointfile,['PHI_' upper(ver_point)],'first_sample',samp_min,'num_samples',samp_max-samp_min);
az_bore=getdata_double_channel_piece(pointfile,['AZ_' upper(ver_point)],'first_sample',samp_min,'num_samples',samp_max-samp_min);
el_bore=getdata_double_channel_piece(pointfile,['EL_' upper(ver_point)],'first_sample',samp_min,'num_samples',samp_max-samp_min);
phi_bore(phi_bore<0)+=360;  % Make phi all positive
ct = getdata_double_channel_piece(timefile,['TIME_' upper(ver_time)],'first_sample',samp_min,'num_samples',samp_max-samp_min);


% Allocate some shit
big_cuts=cell(ndet, 1);
big_raw=zeros(samp_max-samp_min, ndet);
kept_frac=zeros(ndet, 1);

% Loop over detectors and read the data

n_fp=6; %number of spider focal planes
fast_slow=20; %number of fast samples per slow sample
hwp_angles=cell(n_fp,1);

for j=1:ndet,
    % Occasionally inform the human of progress
    if (rem(j,10)==0)
        disp(['Reading data ' num2str(j) ' / ' num2str(ndet)])
    end

    xdet=dets(j,1); rdet=dets(j,2); cdet=dets(j,3);
    det_tag=sprintf('X%1dR%02dC%02d',xdet,rdet,cdet);

    % Read bolometer, flag, and gain dirfiles

    det_raw = getdata_double_channel_piece(datfile,[det_tag '_' upper(ver_clean)],'first_sample',samp_min,'num_samples',samp_max-samp_min);
    flags = getdata_double_channel_piece(flagfile,[det_tag '_' upper(ver_flag_comb)],'first_sample',samp_min,'num_samples',samp_max-samp_min);
    gains = getdata_double_channel_piece(gainfile,[det_tag '_' upper(gain_tag)],'first_sample',samp_min,'num_samples',samp_max-samp_min);
    if isempty(hwp_angles{xdet}),
      disp(['reading HWP angles on ' num2str(xdet)])
      mytag=[ 'X' num2str(xdet) '_ANG_' upper(ver_hwp)];
      disp(mytag)
      istart=floor(samp_min/fast_slow);
      istop=ceil(samp_max/fast_slow);
      disp([istart istop])
      tmp=getdata_double_channel_piece(hwpfile,mytag,'first_sample',istart,'num_samples',istop-istart+1);
      myx=istart:istop;myx=myx'*fast_slow;
      newx=(samp_min:(samp_max-1))';
      hwp_use=interp1(myx,tmp,newx)*pi/180; %yup, we live in radians land
      hwp_angles(xdet)=hwp_use;
    end
    %gains=ones(size(det_raw));
    % Apply gains
    det_raw=det_raw.*gains;

    % Copy detector data into big structure
    big_raw(:,j)=det_raw;
    

    % Convert flags into cut regions, in incomprehensible Jon style
    kept_frac(j)=sum(flags==0)/length(flags);

    tmp=[0;flags;0];
    asdf=diff(tmp==0);
    fdsa=find(asdf);
    assert(rem(numel(fdsa),2)==0);
    cut_regions=[fdsa(1:2:end) fdsa(2:2:end)-1];

    % Deal with jumps, in incomprehensible Jon style
    flags=uint32(flags);
    isjump=0;
    for jj=1:length(jump_bits),
        isjump=isjump+bitget(flags,jump_bits(jj));
    end

    cut_regions(:,end+1)=0;
    for jj=1:size(cut_regions,1),
        flub=sum(isjump(cut_regions(jj,1):cut_regions(jj,2)));
        if flub>0,
            cut_regions(jj,3)=1;
        end
    end
    big_cuts(j)=cut_regions;
end

% Close dirfiles like a responsible adult
close_getdata_file(pointfile);
close_getdata_file(datfile);
close_getdata_file(timefile);
close_getdata_file(flagfile);
close_getdata_file(fjumpfile);
close_getdata_file(jumpfile);
close_getdata_file(gainfile);
close_getdata_file(hwpfile);

% Nuke detectors with low fractions of kept data
ii=kept_frac>min_frac;
big_raw=big_raw(:,ii);
dets=dets(ii,:);
big_cuts=big_cuts(ii);
ndet=size(dets,1);

% Keep track of jumps
njump=zeros(size(big_cuts));
for j=1:length(big_cuts),
    njump(j)=sum(big_cuts{j}(:,3));
end;

% Now deal with detector pointing
offsets = parse_spider_offsets(centroid_file);
ii =~ isnan(offsets(:,4));
offsets = offsets(ii,:);
det_angles=parse_spider_offsets(angle_file);
ii=~isnan(det_angles(:,4));
det_angles=det_angles(ii,:);
if (all(all(det_angles(:,1:3)==offsets(:,1:3)))),
  %in this branch, we can just paste the detector angles in since the detector ordering matches
  offsets(:,end+1)=mean(det_angles(:,4:end),2)*pi/180;
else
  %need to do a better detector alignment, but need to be motivated to write the code first
  assert(1==0)
end

big_ra = zeros(size(big_raw));
big_dec = zeros(size(big_raw));
big_2gamma=zeros(size(big_raw));

for j=1:ndet,
  if (rem(j,10)==0)
    disp(['Detector pointing ' num2str(j) ' / ' num2str(ndet)])
  end
  
  xdet = dets(j,1); rdet=dets(j,2); cdet=dets(j,3);
  ii = find( (offsets(:,1)==xdet)&(offsets(:,2)==rdet)&(offsets(:,3)==cdet));

  daz = offsets(ii,4);
  del = offsets(ii,5);
   
  [ra_detector,dec_detector, sin2psi,cos2psi] = qp_detector_pointing_c(ra_bore, dec_bore, phi_bore, ct, daz,del);
  big_ra(:,j)=ra_detector;
  big_dec(:,j)=dec_detector;
  %xxx - make sure this is correct with HWP, detector definitions etc.
  my2gamma=atan2(sin2psi,cos2psi);
  %big_2gamma(:,j)=+0.5*my2gamma+2*hwp_angles{xdet}+offsets(ii,6); %this line matches Ivan's PA's
  big_2gamma(:,j)=my2gamma+4*hwp_angles{xdet}+2*offsets(ii,6);
  

end

% Subtract low order poly
% xxx do we actually want to do this?
disp('Doing background removal')
mat = legendre_mat(length(det_raw),poly_order);
big_clean = zeros(size(big_raw));
nsamp = size(big_raw,1);
kept_frac = zeros(ndet,1);
wts = zeros(ndet,1);

for j=1:ndet,
  isgood=true(nsamp,1);
  cuts=big_cuts{j};
  for jj=1:size(cuts,1),
    isgood(cuts(jj,1):cuts(jj,2))=false;
  end
  kept_frac(j)=sum(isgood)/length(isgood);
  
  njump=sum(cuts(:,3));
  mymat=mat;
  mymat(:,end+njump)=0;
  cuts=cuts(cuts(:,3)>0,:);
  jumps=round( (cuts(:,1)+cuts(:,2))/2);
  %disp([jumps' sum(isgood)])
  for jj=1:njump,
    mymat(1:jumps(jj),poly_order+jj+1)=0;
    mymat(jumps(jj)+1:end,poly_order+jj+1)=1;
  end
  dat_use=big_raw(isgood,j);
  mymat_small=mymat(isgood,:);

  %find jumps that we don't actually want to fit
  %either because everything before them is cut or everything after them is cut
  tot=sum(mymat_small);
  ii=(tot==0)|(tot==tot(1));
  ii(1)=false;
  mymat=mymat(:,~ii);
  mymat_small=mymat_small(:,~ii);

  lhs=mymat_small'*mymat_small;
  rhs=mymat_small'*dat_use;
  %fitp=inv(lhs)*rhs;
  fitp=invsafe(lhs)*rhs;
  big_clean(:,j)=big_raw(:,j)-mymat*fitp;
  tmp=dat_use-mymat_small*fitp;
  wts(j)=std(tmp);
end

% Save TOD structure
tods=allocate_tod_c();
n=length(det_raw);
set_tod_ndata_c(tods,nsamp);
set_tod_rowcol_c(tods,dets(:,2),dets(:,3)); %'cause why not...
set_tod_pointing_saved(tods,big_ra*pi/180,big_dec*pi/180,big_2gamma); % grownups use radians

set_tod_altaz_c(tods,el_bore*pi/180,az_bore*pi/180);

set_tod_timevec_c(tods,ct);
set_tod_data_saved(tods,big_clean);
set_tod_radec_lims_c(tods);
set_tod_filename(tods,tod_filename);
alloc_tod_cuts_c(tods);

for j=1:ndet,
  mycuts=big_cuts{j};
  for jj=1:size(mycuts,1),
    cuts_extend_c(tods,mycuts(jj,1)-1,mycuts(jj,2),dets(j,2),dets(j,3));
  end
end
