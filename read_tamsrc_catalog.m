function[ra,dec,amps,snr]=read_tamsrc_catalog(fname,varargin)

racol=get_keyval_default('ra',3,varargin{:});
deccol=get_keyval_default('dec',4,varargin{:});
ampcol=get_keyval_default('amps',5,varargin{:});
snrcol=get_keyval_default('snr',7,varargin{:});
mdisp([racol deccol ampcol])

cat=read_text_file_comments(fname);

ra=zeros(numel(cat),1);
dec=zeros(numel(cat),1);
amps=zeros(numel(cat),1);
snr=zeros(numel(cat),1);
for j=1:length(cat),
  tags=strsplit(strtrim(cat{j}),' ',true);  
  %if isempty(tags{1})
  %  tags=tags(2:end);
  %end
  %tags=strsplit(cat{j},sprintf('\t '),true);

  tt=tags{racol};
  if ~isempty(find(tt==':'))
    ra(j)=15*convert_rastr(tt);
  else
    ra(j)=str2num(tt);
  end
  %ra(j)=str2num(tags{racol});


  tt=tags{deccol};
  %if max(find(tt==':'))>0,
  if ~isempty(find(tt==':'))
    dec(j)=convert_rastr(tt);
  else
    dec(j)=str2num(tt);
  end
  %dec(j)=str2num(tags{deccol});
  amps(j)=str2num(tags{ampcol});
  snr(j)=str2num(tags{snrcol});

  %rr=str2num(tags{3})/15;
  %rah=floor(rr);
  %rr=60*(rr-rah);
  %ram=floor(rr);
  %ras=60*(rr-ram);
  %tt=strsplit(tags{2},':',true);
  %amp=str2num(tags{5})/scale_fac;
  %if (str2num(tags{7})>4) %only keep 5-sigma sources
  %  fprintf(fid,'%3d %3d %5.2f %s %s %s %10.3f\n',rah,ram,ras,tt{1},tt{2},tt{3},amp);
  %end
end


function[value]=convert_rastr(tt)
tags=strsplit(strtrim(tt),':',true);
ss=1;
if max(find(tags{1}=='-'))>0
  ss=-1;
end

value=0;
fac=1;
for j=1:length(tags),
  value=value+abs(str2num(tags{j}))/fac;
  fac=fac*60;
end
value=value*ss;