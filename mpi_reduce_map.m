function[map]=mpi_reduce_map(map)
if (strcmp(class(map),'int64'))
  assert(class(map)=='int64');
  mm=skymap2octave(map);
  mm=mpi_allreduce(mm);
  octave2skymap(mm,map);
  return;
end
if isstruct(map)
  if isfield(map,'partition')
    assert(isfield(map,'mapptr'))
    mymap=skymap2octave(map.mapptr);
    if ndims(mymap)==3
      npol=size(mymap,1);
    else
      npol=1;
    end
    big_map=zeros(map.partition.global_nx,map.partition.global_ny);
    ii=map.partition.data_patch;
    jj=map.partition.my_patch;
    if npol==1
      big_map(ii(1):ii(2),ii(3):ii(4))=mymap;
      big_map=mpi_allreduce(big_map);
      mypatch=big_map(jj(1):jj(2),jj(3):jj(4));
    else
      mypatch=zeros(npol,jj(2)-jj(1)+1,jj(4)-jj(3)+1);
      for j=1:npol,
        mm=squeeze(mymap(j,:,:));
        big_map(:,:)=0;
        big_map(ii(1):ii(2),ii(3):ii(4))=mm;
        big_map=mpi_allreduce(big_map);
        mypatch(j,:,:)=big_map(jj(1):jj(2),jj(3):jj(4));
      end
    end
    map.map=mypatch;
  end
end
