function[data]=mapset2tod_octave(mapset,tod,which_tod) 
%pass in a mapset that continas a skymap and correlated noise maps, and get
%the expected tod.  Also need to specify which tod you want.  Need the tod 
%because that contains the pointing info.


map2tod(mapset.skymap.mapptr,tod);

if isfield(mapset,'corrnoise')
    corrnoise=mapset.corrnoise(which_tod);

    if (1)
      if (0)
        for j=1:size(corrnoise.vecs,1),
          vecs2tod(tod,corrnoise.map(:,j),corrnoise.vecs(j,:));
        end
      else
        vecs2tod_blas(tod,corrnoise.map,corrnoise.vecs);
      end
      
    else
      data=get_tod_data(tod);
      data=data+corrnoise2tod(corrnoise);
      push_tod_data(data,tod);
    end
end


if isfield(mapset,'timestreams')
    timestreams=mapset.timestreams(which_tod);
    timestreams2tod_blas(tod,timestreams.map)
end

