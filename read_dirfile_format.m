function[fmt]=read_dirfile_format(todname)

fid=fopen([todname '/format'],'r');
if fid==-1
  error(['missing format on ' todname]);
  return
else
  fclose(fid);
end


mylines=read_text_file_comments([todname '/format']);
nl=length(mylines);
field_names=cell(nl,1);
field_types=cell(nl,1);
vals=cell(nl,1);
for j=1:nl,
  [aa,bb]=strtok(mylines{j});
  field_names(j)={aa};
  [aa,bb]=strtok(bb);
  field_types(j)={aa};
  vals(j)={bb};
end

[aa,bb]=sort(field_names);

fmt.names=field_names(bb);
fmt.types=field_types(bb);
fmt.vals=vals(bb);

