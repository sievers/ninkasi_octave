function[data]=getdata_double_channel(myfile,mychan)
[data,nread]=getdata_double_channel_c(myfile,mychan);
if (nread<numel(data))
  data=data(1:nread);
end

