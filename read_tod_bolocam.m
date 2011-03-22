function[tod]=read_tod_bolocam(fname)
if (iscell(fname)),
  for j=length(fname):-1:1,
    tod(j)=read_tod_bolocam(fname{j});
  end
  tod=int64(tod);
  return
end

[nbolo,nsamp,data,flags,ra,dec,el]=read_raw_bolocam(fname);
tod=allocate_tod_c();
rows=ones(nbolo,1);
cols=(1:nbolo)';
dt=0.1;
tt=dt*(0:(nsamp-1));
set_tod_dt_c(tod,dt);
set_tod_ndata_c(tod,nsamp);
set_tod_timevec_c(tod,tt);
set_tod_rowcol_c(tod,rows,cols);
set_cuts_mustang(tod,flags*1e6,rows,cols);
set_tod_pointing_saved(tod,ra*15*pi/180,dec*pi/180);

set_tod_data_saved(tod,data);
set_tod_radec_lims_c (tod);
tod=int64(tod);