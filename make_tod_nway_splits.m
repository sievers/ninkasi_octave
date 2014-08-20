function[myset]=make_tod_nway_splits(tod_names,varargin)
nway=get_keyval_default('nway',2,varargin{:});
len=get_keyval_default('len',86400,varargin{:});
ct=get_tod_ctimes_from_names(tod_names);

tt=ct/len;tt=tt-floor(tt);
vec=sort([tt;1+tt]);
[a,b]=max(diff(vec));
%disp(['max spacing is at ' num2str([a vec(b) vec(b+1)])]);
thresh=0.5*sum(vec(b:b+1));
myday=ct/len-thresh;
myday=myday-floor(min(myday));
crud=myday-floor(myday);disp(['min/max day fractions are ' num2str([min(crud) max(crud)])])
myset=rem(floor(myday),nway);





