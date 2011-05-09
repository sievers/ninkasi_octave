function[mapset]=tod2mapset_octave(mapset,tod,which_tod)


if isfield(mapset,'skymap')
  tod2map(tod,mapset.skymap.mapptr);
  mapset.skymap.map=skymap2octave(mapset.skymap.mapptr);
end

if isfield(mapset,'corrnoise'),

  if (1)
    if (1)
      %clear_cut_data(tod); %don't use this.
      mapset.corrnoise(which_tod).map=tod2vecs_blas(tod,mapset.corrnoise(which_tod).vecs);
    else      
      for j=1:size(mapset.corrnoise(which_tod).vecs,1),
        mapset.corrnoise(which_tod).map(:,j)=tod2vecs(tod, mapset.corrnoise(which_tod).vecs(j,:));
      end
    end
  else
    mapset.corrnoise(which_tod)=tod2corrnoise(get_tod_data(tod),mapset.corrnoise(which_tod));
  end
  if (0)
    tic
      %clear_cut_data(tod);
      crap=tod2vecs_blas(tod,mapset.corrnoise(which_tod).vecs);
    toc
    
    nn=round(get_tod_ndata(tod)/2);
    crud=mapset.corrnoise(which_tod).map;
    %crap(nn+(1:3),1:3)
    %crud(nn+(1:3),1:3)
    %whos
    disp([ sum(sum(abs(crap-mapset.corrnoise(which_tod).map))) sum(sum(abs(mapset.corrnoise(which_tod).map)))])
  end
end

if isfield(mapset,'timestreams'),
  mapset.timestreams(which_tod).map=tod2timestreams_blas(tod);
end
if isfield(mapset,'srccat')
  if iscell(mapset.srccat),
    for ss=1:numel(mapset.srccat)
      mapset.srccat(ss)={tod2srccat(tod,mapset.srccat{ss})};
    end
  else
    mapset.srccat=tod2srccat(tod,mapset.srccat);
  end
end
