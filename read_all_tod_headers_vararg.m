function[tods,lims,all_lims,isrising]=read_all_tod_headers_vararg(todlist,varargin)


ntod=length(todlist);
tods=int64(zeros(ntod,1));
lims=zeros(ntod,4);
tod_ok=true(size(todlist));
isrising=false(size(todlist));



decimate=get_keyval_default('decimate',0,varargin{:});
pfiles=get_keyval_default('pointing',{},varargin{:});
fix_altaz=get_keyval_default('fix_altaz',true,varargin{:});
tau_files=get_keyval_default('tau_files',[],varargin{:});
isf3db=get_keyval_default('is3db',true,varargin{:});
if iscell(tau_files),
  for j=1:length(tau_files),
    if ~isempty(tau_files{j})
      tmp=load(tau_files{j});
      if (~isf3db)
        tmp(tmp>0)=1.0./(2*pi*tmp(tmp>0));
      end
      %tau_files(j)={load(tau_files{j})};
      tau_files(j)={tmp};

    end
  end
end


assert(~isempty(pfiles));
if iscell(pfiles),
  poffs=cell(size(pfiles));
  for j=1:length(pfiles)
    if ~isempty(pfiles{j})
      poffs(j)=read_pointing_offsets(pfiles{j});
    end
  end
else
  assert(ischar(pfiles));
  poffs=read_pointing_offsets(pfiles);
end
%mdisp('starting tod header read');
for j=1:ntod,
  
  tod=read_tod_header_nopoint_c(todlist{j},decimate);
  if (fix_altaz)
    repair_tod_altaz(tod);
  end
  %mdisp('read raw header');
  
  if iscell(poffs),
    my_poff=poffs{guess_tod_season(todlist{j})};
  else
    my_poff=poffs;
  end 

  %mdisp('setting pointing')
  [mylims,myrising]=set_tod_pointing_100(tod,my_poff);

  %mdisp('read set pointing');
  tods(j)=tod;
  isrising(j)=myrising;
  lims(j,:)=mylims';
end

mdisp('getting to time constant barrier');
mpi_barrier
mdisp('past time constant barrier');

if isempty(tau_files)
  assign_tod_time_constants(tods);
else
  assign_tod_time_constants(tods,tau_files);
end

mpi_barrier
mdisp('past time constant read');


%mdisp('set time constants');
all_lims=lims;
lims=[min(lims(:,1)) max(lims(:,2)) min(lims(:,3)) max(lims(:,4))];

return







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



