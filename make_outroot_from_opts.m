function[froot]=make_outroot_from_opts(froot,myopts)
myopts=set_default_mapping_opts(myopts);



if ~isempty(myopts.noise_scale_facs),
  ind=1:length(myopts.noise_scale_facs);
  ii=ind(myopts.noise_scale_facs==0);
  if ~isempty(ii)
    ii=max(ii);
    froot=[froot '_downweight_' num2str(myopts.bands(ii+1))];
  end
end
if myopts.find_modes_new,
  froot=[froot '_find_new_modes'];
  for j=1:length(myopts.mode_thresh),
    froot=[froot '_' num2str(myopts.mode_thresh(j))];
  end
else
  if myopts.nbadmode>0,
    froot=[froot '_find_' num2str(myopts.nbadmode) '_modes'];
  end
end

if myopts.deconvolve_tau, froot=[froot '_detau'];end;
if myopts.hilton_noise, froot=[froot '_hilton_noise_seed_' num2str(myopts.seed) ];end;
if myopts.do_noise, froot=[froot '_noise'];end;
if myopts.dedark, froot=[froot '_dedark'];end;
if myopts.debutter, froot=[froot '_debutter'];end;
if (myopts.prior_fac>0)
  froot=[froot '_prior_fac_x' num2str(myopts.prior_fac) ];
else
  froot=[froot '_noprior'];
end
if isfield(myopts,'rlcut')
  froot=[froot '_rlcut_' myopts.rlcut];
end
if isfield(myopts,'cut_magic_carpets')
  if myopts.cut_magic_carpets
    froot=[froot '_no_carpets'];
  end
end

if isfield(myopts,'highpass')
  if myopts.highpass>0
    froot=[froot '_highpass_' num2str(myopts.highpass)];
  end
end


if isfield(myopts,'remove_ungainly')
  if myopts.remove_ungainly
    froot=[froot '_gainly'];
  end
end

if isfield(myopts,'remove_common')
  if myopts.remove_common,
    froot=[froot '_uncommon'];
  end
end





froot=[froot '_'];
