function[az_new,tvec_new]=rescale_tod_housekeeping(tod,varargin)

dt=get_keyval_default('dt',1/400,varargin{:});
scan_length=get_keyval_default('scan_length',0,varargin{:});
time_fac=get_keyval_default('time_fac',1,varargin{:});
scan_width=get_keyval_default('scan_width',15,varargin{:});
scan_speed=get_keyval_default('scan_speed',2,varargin{:});
del=get_keyval_default('del',0,varargin{:});
daz=get_keyval_default('daz',0,varargin{:});

[alt,az]=get_tod_altaz(tod);
tvec=get_tod_tvec(tod);

t_cent=median(tvec);
if scan_length==0
  scan_length=tvec(end)-tvec(1);
  scan_length=scan_length*time_fac;
end
tvec_new=((t_cent-0.5*scan_length+dt):dt:(t_cent+0.5*scan_length))';

daz=scan_speed*dt;daz=daz*pi/scan_width;
az_new=(1:length(tvec_new))*daz;az_new=asin(sin(az_new'))*scan_width/(2*pi/2);
mdisp(['scan width is ' num2str(max(az_new)-min(az_new)) ' with step ' num2str(median(abs(diff(az_new))))]);

az_new=az_new*pi/180+median(az)+daz;
free_tod_altaz_c(tod);
free_tod_timevec_c(tod);
free_tod_hwp_c(tod);  

alt_new=1e-7*randn(size(az_new))+median(alt)+del;

set_tod_ndata_c(tod,length(az_new));
set_tod_altaz_c(tod,alt_new,az_new);
set_tod_timevec_c(tod,tvec_new);
set_tod_dt_c(tod,dt);


disp(['initializing actpol pointing.'])
initialize_actpol_pointing(tod);
disp(['initialized.'])




