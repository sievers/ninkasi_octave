function[data]=getdata_double_channel(myfile,mychan)

if ischar(myfile)
  fname=[myfile '/' mychan];
  fid=fopen(fname,'r');
  if fid==-1
    %disp(['error reading ' fname])
    fname2=[myfile '/hk/' mychan];
    fid=fopen(fname2,'r');
    if fid==-1
      error(['cannot find either ' fname ' or ' fname2 ]);
    end
  end
  myprecision='int';
  if strcmp(mychan,'sync_time')
    myprecision='double';
  end
  data=fread(fid,inf,myprecision);
  fclose(fid);
  return
end


if numel(myfile)>1
  data=[];
  for j=1:length(myfile)
    data=[data;getdata_double_channel(myfile(j),mychan)];
  end;
  return;
end

[data,nread]=getdata_double_channel_c(myfile,mychan);
if (nread<numel(data))
  data=data(1:nread);
end

