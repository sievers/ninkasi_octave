function[mapset]=apply_precon(mapset,precon)
global tods

do_cutvecs=false;
if ~isempty(tods)
  do_cutvecs=are_cutvecs_fitparams(tods(1));
  for j=2:numel(tods),
    assert(do_cutvecs==are_cutvecs_fitparams(tods(j)));
  end
end

if isfield(mapset,'skymap')
  if isfield(precon,'naess_precon')    
    mapset.skymap=apply_sigurd_precon(mapset.skymap,precon.naess_precon);
  else
    if isfield(precon,'skymap'),  %apply a Jacoby preconditioner to the sky
      if (~is_map_polarized(mapset.skymap.mapptr))
        mapset.skymap.map(precon.skymap.map>0)=mapset.skymap.map(precon.skymap.map>0)./precon.skymap.map(precon.skymap.map>0);
      else
        if isfield(mapset.skymap,'partition')
          %have to make a map here because the C version of the preconditioner
          %will cover a different patch of sky than the C version of the regular map
          %(but should be the same region as the octave version of the regular map)

          tmp=make_map_copy(precon.skymap.mapptr);
          set_map_polstate_c(tmp,get_map_polstate_c(mapset.skymap.mapptr));
          octave2skymap(mapset.skymap.map,tmp);
          apply_pol_precon_c(tmp,precon.skymap.mapptr);

          mapset.skymap.map=skymap2octave(tmp);
          destroy_map(tmp);

          
          
        else
          octave2skymap(mapset.skymap);
          apply_pol_precon_c(mapset.skymap.mapptr,precon.skymap.mapptr);
          mapset.skymap.map=skymap2octave(mapset.skymap.mapptr);
        end
      end
      
      %mapset.skymap.map=mapset.skymap.map.*precon.skymap.map;
    end
  end
end

if isfield(mapset,'jumps'),
  if isfield(precon,'jumps'),
    if iscell(mapset.jumps),
      for j=1:length(mapset.jumps),
        mapset.jumps{j}=mapset.jumps{j}./precon.jumps{j};
      end
    else
      mapset.jumps=mapset.jumps./precon.jumps;
    end
  end
end


if isfield(mapset,'ground')
  if isfield(precon,'ground')
    
    octave2skymap(precon.ground.ground,precon.ground.groundptr);
    octave2skymap(mapset.ground.ground,mapset.ground.groundptr);
    apply_pol_precon_c(mapset.ground.groundptr,precon.ground.groundptr);
    mapset.ground.ground=skymap2octave(mapset.ground.groundptr);
  end
end



%since we often cut windows, then use a preconditioner for the cut section.
if isfield(mapset,'cutvecs')
  if isfield(precon,'cutvecs')
    %mdisp('applying cutvecs preconditioner')
    if iscell(mapset.cutvecs)
      for j=1:numel(mapset.cutvecs),
        if (do_cutvecs)
          %mdisp('calling fitparams precon')
          mapset.cutvecs(j)={apply_cutvecs_fitparams_precon(tods(j),mapset.cutvecs{j})};
        else
          mapset.cutvecs(j)={mapset.cutvecs{j}./precon.cutvecs{j}};
        end
      end
    else
      if (do_cutvecs)
        %mdisp('calling fitparams precon')
        mapset.cutvecs=apply_cutvecs_fitparams_precon(tods,mapset.cutvecs);
      else
        mapset.cutvecs=mapset.cutvecs./precon.cutvecs;
      end
    end
  end
end

                      

return

if isfield(precon,'corrnoise')
  for j=1:length(precon.corrnoise),
    if (isfield(precon(j).corrnoise,'prior'))
      mapset.corrnoise(j).map= fit_corrnoise_wpriors(precon.corrnoise(j).vecs,mapset.corrnoise(j).map,precon.corrnoise(j).prior);
    else
      mapset.corrnoise(j).map= fit_corrnoise_wpriors(precon.corrnoise(j).vecs,mapset.corrnoise(j).map);
    end
  end
end

      