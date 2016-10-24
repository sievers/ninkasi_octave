function[mapset]=create_prefiltered_map(tods,mapset,mapset_in,myopts)
mapset=clear_mapset(mapset);
do_actpol_pointing=get_struct_mem(myopts,'do_actpol_pointing',false);
free_2gamma=get_struct_mem(myopts,'free_2gamma',1);


for j=1:length(tods),
  mytod=tods(j);
  create_initial_mapset_opts(mytod,mapset,mapset_in,myopts);
  

  if (do_actpol_pointing)

    precalc_actpol_pointing_exact(mytod,1);
    if isfield(mapset.skymap,'mapptr')
      convert_saved_pointing_to_pixellization(mytod,mapset.skymap.mapptr)
    end
    free_tod_pointing_saved(mytod,free_2gamma);

    if (~does_tod_have_twogamma_fit(mytod))
      precalc_actpol_pointing_exact(mytod,2);
      set_tod_twogamma_fit(mytod,'npoly_2gamma',3);
      free_tod_pointing_saved(mytod,free_2gamma);
    end

  else
    %for non-actpol stuff, make sure this gets set up sensibly.
    assert(1==0);
  end
  
  mapset=tod2mapset_octave(mapset,mytod,j);
  free_tod_storage(mytod);
  if (do_actpol_pointing) %we cached the pointing earlier, now we need to get rid of it.
    free_tod_pointing_saved(mytod,free_2gamma); %should have already been freed, but just in case...
    free_saved_pixellization(mytod);
  end
end

if isfield(mapset,'skymap')
  %do the mpi reduction
  mapset.skymap.map=skymap2octave(mapset.skymap.mapptr);
  mapset.skymap.map=mpi_allreduce(mapset.skymap.map);
  octave2skymap(mapset.skymap);
end
