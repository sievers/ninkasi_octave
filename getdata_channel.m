function[value]=getdata_channel(myf,fmt,fname)
ii=_find_dirfile_tag(fmt,fname);
%disp(['type is ' fmt.types{ii}]);
%disp(['val  is ' fmt.vals{ii}]);
if ii==-1
  value=[];
  return
end


if strcmp(fmt.types{ii},'RAW')
  value=getdata_double_channel(myf,fname);
  return
end

if strcmp(fmt.types{ii},'LINCOM')
  value=0;
  tags=strsplit(fmt.vals{ii},' ',true);
  nterm=str2num(tags{1});
  tags=tags(2:end);
  for j=1:nterm,
    value=value+str2num(tags{3*j})+str2num(tags{3*j-1})*getdata_channel(myf,fmt,tags{3*j-2});
  end
  
  return
end



function[ind]=_find_dirfile_tag(fmt,fname)
vec=fmt.names;
vec(end+1)={fname};
[aa,bb]=sort(vec);
myind=find(bb==length(bb));
%disp(['possible targets are ' aa{myind-1} ' ' aa{myind+1} ]);
if strcmp(fname,aa{myind-1})==1
  ind=myind-1;
  return
end
if strcmp(fname,aa{myind+1})==1
  ind=myind+1;
  return
end
error(['key ' fname ' not found in format file.']);
ind=-1;
return
