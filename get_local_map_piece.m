function[my_patch]=get_local_map_piece(map)
assert(isstruct(map));
assert(isfield(map,'partition'));
assert(isfield(map,'map'));

if ndims(map.map)==3
  npol=size(map.map,1);
else
  npol=1;
end


nproc=mpi_comm_size;
myid=mpi_comm_rank+1;

big_map=zeros(map.partition.global_nx,map.partition.global_ny);
jj=map.partition.my_patch;
ii=map.partition.data_patch;

if npol==1,
  big_map(jj(1):jj(2),jj(3):jj(4))=map.map;
  big_map=mpi_allreduce(big_map);
  my_patch=big_map(ii(1):ii(2),ii(3):ii(4));
else
  mypatch=zeros(npol,ii(2)-ii(1)+1,ii(4)-ii(3)+1);
  for j=1:npol,
    big_map(:,:)=0;
    big_map(jj(1):jj(2),jj(3):jj(4))=squeeze(map.map(j,:,:));
    big_map=mpi_allreduce(big_map);
    my_patch(j,:,:)=big_map(ii(1):ii(2),ii(3):ii(4));
  end
end

    


  