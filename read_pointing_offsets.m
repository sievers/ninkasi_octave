function[poff]=read_pointing_offsets(fname,varargin)

flip_az=get_keyval_default('flip_az',false,varargin{:});
flip_alt=get_keyval_default('flip_alt',false,varargin{:});
az_shift=get_keyval_default('az_shift',0,varargin{:});
alt_shift=get_keyval_default('alt_shift',0,varargin{:});



fid=fopen(fname,'r');
if fid==-1
  warning(['Unable to read file ' fname ' in read_pointing_offsets.']);
  poff=[];
  return
end



ll=fgetl(fid);
args=line_to_argv(ll);
nrow=str2num(args{3});
ncol=str2num(args{6});
ll=fgetl(fid);
args=line_to_argv(ll);
fittype=str2num(args{3});
if fittype==100
  lines=read_text_file_comments(fid,'#');
  fclose(fid);
  poff=parse_pointing_offsets_type_100(lines);
  poff.nrow=nrow;
  poff.ncol=ncol;
  return
end


vals=fscanf(fid,'%f',[2 inf]);
fclose(fid);

assert(length(vals)==nrow*ncol+1);
assert((fittype==0)|(fittype==10));
alt_off=reshape(vals(1,2:end)',[ncol nrow])';
az_off=reshape(vals(2,2:end)',[ncol nrow])';

if ~isempty(findstr(fname,'rising'))
  if az_shift~=0
    disp(['shifting az by ' num2str(az_shift) ' on ' fname]);
    az_off=az_off+az_shift;
  end
  
  if alt_shift~=0
    disp(['shifting alt by ' num2str(alt_shift) ' on ' fname]);
    alt_off=alt_off+alt_shift;
  end
end




if (flip_az|flip_alt)
  vec=reshape(az_off,[numel(az_off) 1]);
  vec=vec(abs(vec)<10);
  az_cent=mean([min(vec) max(vec)]);
  %az_cent=median(vec(abs(vec)<10));
  vec=reshape(alt_off,[numel(alt_off) 1]);
  vec=vec(abs(vec)<10);
  alt_cent=mean([min(vec) max(vec)]); 
  %alt_cent=median(vec(abs(vec)<10))
  %disp(['centers are at ' num2str([az_cent alt_cent])])
  ind=abs(alt_off)<10;
  if flip_az,
    az_off(ind)=2*az_cent-az_off(ind);
  end
  if (flip_alt)
    alt_off(ind)=2*alt_cent-alt_off(ind);
  end
end




poff.dalt=alt_off;
poff.daz=az_off;
poff.fittype=fittype;
poff.fitparams=vals(:,1);




function[myargv]=line_to_argv(ll)
myargv={};
ii=0;
while(1)
  [tag,ll]=strtok(ll);
  ii=ii+1;
  myargv(ii)={tag};
  if length(ll)==0
    return;
  end
end


