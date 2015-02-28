function[taus]=get_pertod_taus(tod,varargin)
dirroot=get_keyval_default('dirroot','/home/r/rbond/egrace/depot/tau_ar1_v2/',varargin{:});
tailtag=get_keyval_default('tailtag','.ar1.tau',varargin{:});
skip_ct2=get_keyval_default('skip_ct2',false,varargin{:});
nrow=get_keyval_default('nrow',32,varargin{:});
ncol=get_keyval_default('ncol',32,varargin{:});
todname=get_tod_name(tod);
if iscell(todname)
  todname=todname{1};
  end

[ct1,ct2]=get_tod_ctimes_from_names(todname);
dr2=sprintf('%d',ct1);dr2=dr2(1:5);
if skip_ct2
  tauname=[dirroot '/' dr2 '/' sprintf('%d',ct1) tailtag];
else
  tauname=[dirroot '/' dr2 '/' sprintf('%d',ct1) '.' sprintf('%d',ct2) tailtag];
end



tau_raw=load(tauname);

taus=zeros(nrow+1,ncol+1);

rr=floor(tau_raw(:,1)/ncol);
cc=tau_raw(:,1)-rr*ncol;

for j=1:length(tau_raw),
  taus(rr(j)+1,cc(j)+1)=tau_raw(j,2);
end



