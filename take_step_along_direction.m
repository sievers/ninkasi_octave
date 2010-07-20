function[new_x,a]=take_step_along_direction(dir,x,b,tods)
%function[value]=take_step_along_direction(dir,x,b).
%find a such that chi^2 is minimized for x+a*dir

Adir=mapset2mapset_corrnoise_octave(tods,dir);

lhs=mapsetdotmapset(Adir,dir);

rhs=mapsetdotmapset(dir,b)-mapsetdotmapset(Adir,x);


lhs_prior=0;
rhs_prior=0;
for j=1:length(tods)
  if isfield(dir.corrnoise(j),'prior'),
    gam_q=apply_corrnoise_priors(x.corrnoise(j).vecs, dir.corrnoise(j).map,x.corrnoise(j).prior);
    lhs_prior=lhs_prior+sum(sum(dir.corrnoise(j).map.*gam_q));
    rhs_prior=rhs_prior+sum(sum(x.corrnoise(j).map.*gam_q));
  end
end

lhs_prior=mpi_allreduce(lhs_prior);
rhs_prior=mpi_allreduce(rhs_prior);

a=(rhs-rhs_prior)/(lhs+lhs_prior);
new_x=mapset_axpy(x,dir,a);


