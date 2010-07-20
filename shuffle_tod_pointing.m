function[value]=shuffle_tod_pointing(tods)
for j=1:length(tods),
  mylims(j,:)=get_tod_radec_lims_c(tods(j));
  altaz=get_tod_altaz_lims_c(tods(j));
  myaz(j,1)=mean(altaz(3:4));
  local_lims=[min(mylims(:,1)) max(mylims(:,2)) min(mylims(:,3)) max(mylims(:,4))];
end

global_lims(1)=mpi_allreduce(local_lims(1),'min');
global_lims(2)=mpi_allreduce(local_lims(2),'max');
global_lims(3)=mpi_allreduce(local_lims(3),'min');
global_lims(4)=mpi_allreduce(local_lims(4),'max');

rising=myaz<pi;

shuffle_tod_set(tods(rising),mylims(rising,:),global_lims);
shuffle_tod_set(tods(~rising),mylims(~rising,:),global_lims);


return


function[value]=shuffle_tod_set(tods,mylims,global_lims)


ntods=mpi_concatenate(numel(tods));
myrank=mpi_comm_rank+1;
myind=(1+sum(ntods(1:myrank-1))):sum(ntods(1:myrank));
nproc=mpi_comm_size;

all_lims=mpi_concatenate(mylims);

ntod=sum(ntods);
if (myrank==1)
  to_shuffle=randperm(ntod);
else
  to_shuffle='';
end
to_shuffle=mpi_bcast_array(to_shuffle,0);
my_new_lims=all_lims(to_shuffle(myind),:);
assert(size(my_new_lims,1)==length(tods));

%mdisp(['global lims are ' num2str(global_lims)])
%for jj=1:nproc,
%  nn=mpi_allreduce(1);
%  if (myrank==jj)
for j=1:length(tods),
  [ra_shift,dec_shift]=get_inbounds_shifts(mylims(j,:),my_new_lims(j,:),global_lims);  
  %ra_shift=1e-6*ra_shift;
  %dec_shift=0;
  shift_tod_pointing_c(tods(j),ra_shift,dec_shift*0);
  disp(['shifting tod w/ lims ' num2str(mylims(j,:)) ' from ' num2str(my_new_lims(j,:)) ' by ' num2str([ra_shift dec_shift])]);
end
%  end
%  pause(0.1);
%end




return



function[ra_shift,dec_shift]=get_inbounds_shifts(mylims,my_new_lims,global_lims)
ra_shift=mean(my_new_lims(1:2))-mean(mylims(1:2));

if (ra_shift+mylims(2)>global_lims(2))
  ra_shift=global_lims(2)-mylims(2);
end

if (ra_shift+mylims(1)<global_lims(1))
  ra_shift=global_lims(1)-mylims(1);
end

dec_shift=mean(my_new_lims(3:4))-mean(mylims(3:4));

if (dec_shift+mylims(4)>global_lims(4))
  dec_shift=global_lims(4)-mylims(4);
end

if (dec_shift+mylims(3)<global_lims(3))
  dec_shift=global_lims(3)-mylims(3);
end








