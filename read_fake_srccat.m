function[ra,dec,amps]=read_fake_srccat(fname)
lines=read_text_file_comments(fname);
n=numel(lines)
ra=zeros(n,1);
dec=zeros(n,1);
amps=zeros(n,1);
for j=1:n,
  tags=strsplit(strtrim(lines{j}),' ',true);
  nn=length(tags);
  if nn==3,
    ra(j)=str2num(tags{1});
    dec(j)=str2num(tags{2});
    amps(j)=str2num(tags{3});
  else
    v=zeros(nn,1);
    for k=1:nn,
      v(k)=abs(str2num(tags{k}));
    end
    mysign=1;
    if ~isempty(find(tags{4}=='-'))
      mysign=-1;
    end
    ra(j)=15*(v(1)+v(2)/60+v(3)/3600);
    dec(j)=mysign*(v(4)+v(5)/60+v(6)/3600);
    amps(j)=v(7);
  end
end
