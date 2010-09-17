function[data,vecs,fitp]=get_tod_data_cleaned(tod,varargin)

vecs=get_keyval_default('vecs',[],varargin{:});
do_cm=get_keyval_default('cm',true,varargin{:});
imin=get_keyval_default('imin',1,varargin{:});
imax=get_keyval_default('imax',get_tod_ndata(tod),varargin{:});
order=get_keyval_default('order',-1,varargin{:});
freqs=get_keyval_default('freqs',[],varargin{:});
do_detrend=get_keyval_default('do_detrend',true,varargin{:});

%be able to specify incoming data different from the TOD if you wish.
%say, when modelling the correlated timestreams.
data=get_keyval_default('data',[],varargin{:});
if isempty(data)
  if do_detrend    
    data_org=get_tod_data(tod);
    gapfill_data_c(tod);
    array_detrend(tod);
    gapfill_data_c(tod);
    window_data_c(tod);
    data=get_tod_data(tod);
    push_tod_data(data_org,tod);
    clear data_org
  else
    data=get_tod_data(tod);
  end
else
  disp('data coming in from outside')
end

tvec=get_tod_tvec(tod);
tt=tvec-mean(tvec);
tt=tt/max(abs(tt));
for j=0:order,
  vecs=[vecs tt.^j];
end
%for j=1:length(freqs),
  %vecs=[vecs cos(2*pi*freqs*tvec) sin(2*pi*freqs*tvec)];
  vecs=[vecs cos(2*pi*tvec*freqs) sin(2*pi*tvec*freqs)];
%end
if(do_cm)
  vecs=[vecs get_cleaned_cm(data)];
end
rhs=vecs(imin:imax,:)'*data(imin:imax,:);
fitp=inv(vecs(imin:imax,:)'*vecs(imin:imax,:))*rhs;
data=data-vecs*fitp;


if imin>1,
  data(1:imin,:)=repmat(data(imin+1,:),[imin 1]);
end


nn=get_tod_ndata(tod);
if imax<nn
  data(imax:nn,:)=repmat(data(imax-1,:),[(nn-imax)+1 1]);
end


return



function[cm]=get_cleaned_cm(data)
cm=median(data,2);
amps=(cm'*data)/(cm'*cm);
for j=1:size(data,2),
  data(:,j)=data(:,j)/amps(j);
end
cm=median(data,2);