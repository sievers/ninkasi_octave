function[data]=mapset2tod_octave(mapset,tod,which_tod) 
%pass in a mapset that continas a skymap and correlated noise maps, and get
%the expected tod.  Also need to specify which tod you want.  Need the tod 
%because that contains the pointing info.



[myptr,mytype]=get_generic_tod_pointer(tod);

if isfield(mapset,'skymap')
  if (mytype==0)
    map2tod(mapset.skymap.mapptr,tod);
  end
  if (mytype==1)
    map_to_vis(myptr,mapset.skymap.mapptr);
  end
end


if (mytype==0)
  if isfield(mapset,'cutvecs')
    if iscell(mapset.cutvecs),    
      cutvec2tod_c(tod,mapset.cutvecs{which_tod});
    else
      cutvec2tod_c(tod,mapset.cutvecs);
    end
  end
end

if (mytype==0)
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
end

if mytype==0
  if isfield(mapset,'timestreams')
    timestreams=mapset.timestreams(which_tod);
    timestreams2tod_blas(tod,timestreams.map)
  end
end




if isfield(mapset,'srccat')
  if (mytype==0)   
    if iscell(mapset.srccat),
      for ss=1:numel(mapset.srccat)
        add_srccat2tod(tod,mapset.srccat{j});
      end
    else
      add_srccat2tod(tod,mapset.srccat);
    end
  else
    warning(['sources being ignored on tod type ' num2str(mytype)])
  end
end


