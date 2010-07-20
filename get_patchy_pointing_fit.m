function[altvec,azvec,ra_mat,dec_mat,ra_clock,dec_clock]=get_patchy_pointing_fit(tod,delta,myord,ra_wrap)
altaz_lims=get_tod_alldet_altaz_lims_c(tod);
if ~exist('delta')
  delta=0.002;
end
if ~exist('myord')
  myord=2;
end
if ~exist('ra_wrap')
  ra_wrap=pi;
end
mdisp(['ra_wrap is ' num2str(ra_wrap)]);
az1=altaz_lims(3)-0.01*(altaz_lims(4)-altaz_lims(3));azvec=az1:delta:(altaz_lims(4)+delta);
alt1=altaz_lims(1)-0.01*(altaz_lims(2)-altaz_lims(1));altvec=alt1:delta:(altaz_lims(2)+delta);

altmat=repmat(altvec',[1 length(azvec)]);
azmat=repmat(azvec,[length(altvec) 1]);
altmat_as_vec=reshape(altmat,[numel(altmat) 1]);
azmat_as_vec=reshape(azmat,[numel(azmat) 1]);

ctime_ref=get_tod_ctime_c(tod)+0.5*(get_tod_dt(tod)*get_tod_ndata(tod));

ctimes=zeros(size(azmat_as_vec))+ctime_ref;
[ra_as_vec,dec_as_vec]=get_radec_from_altaz_ctime_c (altmat_as_vec,azmat_as_vec,ctimes);

ra_as_vec=unwrap_ra(ra_as_vec,ra_wrap);

ra_mat=reshape(ra_as_vec,size(altmat));
dec_mat=reshape(dec_as_vec,size(altmat));

ndata=get_tod_ndata(tod);
dt=get_tod_dt(tod);


if (1)
  az_cent=median(azvec);
  alt_cent=median(altvec);
  ctime1=get_tod_ctime_c(tod);
  ctime_vec=ctime1:(ctime1+dt*ndata);  %evaluate every second
  ctime_vec=[ctime_ref ctime_vec];
  [ra_as_vec2,dec_as_vec2]=get_radec_from_altaz_ctime_c (alt_cent+0*ctime_vec,az_cent+0*ctime_vec,ctime_vec);
  ra_as_vec2=unwrap_ra(ra_as_vec2,ra_wrap);
  ra_clock=polyfit(ctime_vec-ctime1,ra_as_vec2-ra_as_vec2(1),myord);
  dec_clock=polyfit(ctime_vec-ctime1,dec_as_vec2-dec_as_vec2(1),myord);
else
  [ra_as_vec2,dec_as_vec2]=get_radec_from_altaz_ctime_c (altmat_as_vec,azmat_as_vec,ctimes+dt*ndata);
  ra_clock=median(ra_as_vec2-ra_as_vec)/(dt*ndata);
  dec_clock=median(dec_as_vec2-dec_as_vec)/(dt*ndata);
end


function[ra_new]=unwrap_ra(ra,ra_wrap)
ind=(ra>ra_wrap);

ra_new=ra;
ra_new(ind)=ra_new(ind)-2*pi;

return