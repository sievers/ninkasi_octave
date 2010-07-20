function[vec,left_val,right_val]=array_detrend(tod,varargin)
%subtract a median and a median slope, where the median slope is fixed for all detectors.
assert(class(tod)=='int64');
window_len=get_keyval_default('window_len',20,varargin{:});

dt=get_tod_dt(tod);
nsamp=ceil(window_len/dt);
data=get_tod_data(tod);
%mdisp(['nsamp is ' num2str(nsamp)]);

mymed=omp_median(data,1);
data=data-repmat(mymed,[size(data,1) 1]);
left_val=median(omp_median(data(1:nsamp,:),1));
right_val=median(omp_median(data(end-nsamp+1:end,:),1));
%mdisp(['left/right vals are ' num2str([left_val right_val])]);
nn=size(data,1);
ind=(1:nn)'-round(nsamp/2);ind=ind/(nn-nsamp);ind=ind-mean(ind);
vec=ind*(right_val-left_val);
data=data-repmat(ind*(right_val-left_val),[1 size(data,2)]);
push_tod_data(data,tod);
