function[tod]=read_tod_dooley(froot,varargin)

which_pol=get_keyval_default('pol','I',varargin{:});



if iscell(froot)

  for j=1:length(fname)
    %tod(j)=read_tod_mustang(fname{j});
    if j>1
      tod=[tod read_tod_dooley(froot{j},varargin{:})];
    else
      tod=read_tod_dooley(froot{j},varargin{:});
    end
  end
  return
end


[ex,ey,ra,dec,pa,tt]=read_dooley(froot);
tod=allocate_tod_c();
rows=1:size(ex,2);
cols=ones(size(rows));
dt=median(diff(tt));


set_tod_dt_c(tod,dt);
set_tod_ndata_c(tod,length(tt));
set_tod_timevec_c(tod,tt);
set_tod_rowcol_c(tod,rows,cols);
set_tod_pointing_saved (tod,ra,dec);

alloc_tod_cuts_c(tod);

switch(lower(which_pol))
  case{'i'}
   data=(ex+ey)/4;
 otherwise
  error(['unrecognized polarization ' which_pol ' requested.']);
end
set_tod_data_saved(tod,data);
set_tod_radec_lims_c (tod);
tod=int64(tod);
   