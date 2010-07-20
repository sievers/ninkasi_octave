function[myopts]=set_default_mapping_opts(myopts)
if ~exist('myopts')
  myopts=[];
end

myopts=set_struct_mem_default(myopts,'add_noise',false);
myopts=set_struct_mem_default(myopts,'keep_data',false);
myopts=set_struct_mem_default(myopts,'scale_factor',1);
myopts=set_struct_mem_default(myopts,'detrend',false);
myopts=set_struct_mem_default(myopts,'detrend_array',false);
myopts=set_struct_mem_default(myopts,'debutter',false);
myopts=set_struct_mem_default(myopts,'highpass',0);
myopts=set_struct_mem_default(myopts,'check_for_nans',true);
myopts=set_struct_mem_default(myopts,'deconvolve_tau',false);
myopts=set_struct_mem_default(myopts,'reverse_time',false);
myopts=set_struct_mem_default(myopts,'hilton_noise',false);
myopts=set_struct_mem_default(myopts,'hilton_nmode',24);
myopts=set_struct_mem_default(myopts,'dedark',false);


myopts=set_struct_mem_default(myopts,'do_noise',false);
myopts=set_struct_mem_default(myopts,'bands',[-1 3 300]);
myopts=set_struct_mem_default(myopts,'rots',[0 0]);
myopts=set_struct_mem_default(myopts,'noise_types',[1 0]);
myopts=set_struct_mem_default(myopts,'noise_scale_facs',[]);
myopts=set_struct_mem_default(myopts,'do_calib',false);
myopts=set_struct_mem_default(myopts,'remove_common',false);

myopts=set_struct_mem_default(myopts,'add_input_map',false);
myopts=set_struct_mem_default(myopts,'do_input_mapset',false);

myopts=set_struct_mem_default(myopts,'find_modes',false);
myopts=set_struct_mem_default(myopts,'find_modes_new',false);
myopts=set_struct_mem_default(myopts,'write_cleaned_data',false);

myopts=set_struct_mem_default(myopts,'nu1',0.01);
myopts=set_struct_mem_default(myopts,'nu2',4);
myopts=set_struct_mem_default(myopts,'nu3',30);
myopts=set_struct_mem_default(myopts,'nbadmode',3);

myopts=set_struct_mem_default(myopts,'seed',0);
myopts=set_struct_mem_default(myopts,'prior_fac',0);
myopts=set_struct_mem_default(myopts,'gaussian_noise',false);
myopts=set_struct_mem_default(myopts,'monitor_tods',false);
signal_only=get_struct_mem(myopts,'signal_only',false);







return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[myopts]=set_struct_mem_default(myopts,mystr,myval)
if ~isfield(myopts,mystr)
  eval(['myopts.' mystr ' = myval;']);
end
