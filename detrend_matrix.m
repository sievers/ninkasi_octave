function[data,vec,left_val,right_val]=detrend_matrix(data,varargin)
%subtract a median and a median slope, where the median slope is fixed for all detectors.
nsamp=get_keyval_default('npad',1,varargin{:});

mymed=omp_median(data,1);
data=data-repmat(mymed,[size(data,1) 1]);
left_val=median(omp_median(data(1:nsamp,:),1));
right_val=median(omp_median(data(end-nsamp+1:end,:),1));
nn=size(data,1);
ind=(1:nn)'-round(nsamp/2);ind=ind/(nn-nsamp);ind=ind-mean(ind);
vec=ind*(right_val-left_val);
data=data-repmat(ind*(right_val-left_val),[1 size(data,2)]);
