function[tod_names]=remove_ungainlies(tod_names,varargin)

if 1
  %do varargin stuff here.
  fname=get_keyval_default('fname','/home/sievers/act/blacklists/ungainly_ar1_2008_091102.txt',varargin{:});
  %/home/sievers/act/blacklists/ungainly_ar2_2008_091130.txt
  n1=length(tod_names);
  tod_names=remove_tods_from_blacklist(tod_names,fname);
  n2=length(tod_names);
  mdisp(['removed ' num2str(n1-n2) ' tods in ungainly cuts from file ' fname]);
  
  return
  
end

ungainly_name='ungainly_ar1_2008_091102.txt';
ungainly_dir='/home/sievers/act/blacklists/';

n1=length(tod_names);
tod_names=remove_tods_from_blacklist(tod_names,[ungainly_dir ungainly_name]);
n2=length(tod_names);
mdisp(['removed ' num2str(n1-n2) ' tods in ungainly cuts.']);
