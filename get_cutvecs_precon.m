function[cutvecs]=get_cutvecs_precon(tods,do_symmetric)
if ~exist('do_symmetric')
  do_symmetric=false;
end

if numel(tods)>1,

  cutvecs=cell(size(tods));
  for j=1:numel(tods),
    cutvecs{j}=get_cutvecs_precon(tods(j),do_symmetric);
  end
  return
end


allocate_tod_storage(tods);
assign_tod_value(tods,1);
window_data(tods);
cutvecs=tod2cutvec_c(tods);
cutvecs(cutvecs==0)=1;
if (do_symmetric)
  cutvecs=cutvecs.^2;
end

free_tod_storage(tods);
