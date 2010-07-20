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

  %dr='/cita/d/raid-sievers3/sievers/act/data/hilton_sims/';
  %dr='/cita/d/raid-cita/sievers/act/hilton_sims/';
  dr='/mnt/raid-cita/sievers/act/hilton_sims/';
  todlist={
      [dr '1221390748.1221390761.ar1/']
[dr '1221993788.1221993903.ar1/']
[dr '1222252186.1222252283.ar1/']
[dr '1223630848.1223630964.ar1/']
[dr '1223700209.1234326370.ar1/']
[dr '1223717251.1234327922.ar1/']
[dr '1223975543.1223975689.ar1/']
[dr '1224130828.1224130940.ar1/']
[dr '1224217241.1224217412.ar1/']
[dr '1224475437.1224475573.ar1/']
[dr '1224578536.1224578653.ar1/']
[dr '1224923283.1224923419.ar1/']
[dr '1225009631.1225009755.ar1/']
[dr '1225078521.1225078640.ar1/']
[dr '1225267927.1225268037.ar1/']
[dr '1225354329.1225354462.ar1/']
[dr '1225699022.1234639633.ar1/']
[dr '1225768004.1234649611.ar1/']
[dr '1226112714.1226112873.ar1/']
[dr '1226285106.1226357005.ar1/']
[dr '1226302013.1226348217.ar1/']
[dr '1226560306.1226619322.ar1/']
[dr '1226646707.1226646820.ar1/']
[dr '1226991401.1226991516.ar1/']
[dr '1227491393.1227491515.ar1/']
[dr '1227939086.1227939191.ar1/']
[dr '1228025486.1228025589.ar1/']
[dr '1228283781.1228283887.ar1/']
[dr '1228783817.1228783833.ar1/']
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



