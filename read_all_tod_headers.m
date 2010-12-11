function[tods,lims,all_lims,isrising]=read_all_tod_headers(todlist,decimate,varargin)

if ~exist('todlist')
  todlist='';
end

if isnumeric(todlist)
  if length(todlist)>0
    ntod=todlist;
    todlist='';
  end
end


if isempty(todlist)
  dr='/mnt/raid-cita/sievers/act/hilton_sims/';
  todlist={
      [dr '1221390748.1221390761.ar1/']
          };
  
  
end

if (exist('ntod')),
  todlist=todlist(1:ntod);
end


ntod=length(todlist);
tods=int64(zeros(ntod,1));
lims=zeros(ntod,4);
tod_ok=true(size(todlist));
isrising=false(size(todlist));

for j=1:ntod,
  if exist('decimate')
    [tod,lim,am_i_rising]=read_tod_header(todlist{j},'',decimate,varargin{:});
  else
    [tod,lim,am_i_rising]=read_tod_header(todlist{j},varargin{:});
  end

  isrising(j)=am_i_rising;
  tods(j)=tod;
  if isempty(lim)
    tod_ok(j)=false;
  else
    lims(j,:)=lim;
  end
end
assign_tod_time_constants(tods);

tods=tods(tod_ok);
lims=lims(tod_ok,:);


all_lims=lims;
lims=[min(lims(:,1)) max(lims(:,2)) min(lims(:,3)) max(lims(:,4))];



