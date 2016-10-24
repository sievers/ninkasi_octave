function[left,right]=find_scans_cbass(tod,varargin)

diffthresh=get_keyval_default('diffthresh',0.03*pi/180,varargin{:}); %threshold in daz/d_sample at which we start turning around
scanthresh=get_keyval_default('scanthresh',5000,varargin{:});  %minimum scan length, samples


[alt,az]=get_tod_altaz(tod);
az=unwrap(az)-2*pi;

azdiff=diff(az);

[ibot1,ibot2]=get_scans(azdiff,diffthresh,scanthresh);
[itop1,itop2]=get_scans(-1*azdiff,diffthresh,scanthresh);

left=[ibot1 ibot2];
right=[itop1 itop2];

function[ibot1,ibot2]=get_scans(azdiff,diffthresh,scanthresh)

ibot1=find(diff(azdiff<-diffthresh)==1);
ibot2=find(diff(azdiff<-diffthresh)==-1);

if ibot1(1)>ibot2(1)
  ibot2=ibot2(2:end);
  assert(ibot2(1)>ibot1(1));
end
if numel(ibot1)>numel(ibot2)
  ibot1=ibot1(1:end-1);
  assert(numel(ibot1)==numel(ibot2));
end

scan_len=ibot2-ibot1;
ii=scan_len>scanthresh;
ibot1=ibot1(ii);
ibot2=ibot2(ii);
