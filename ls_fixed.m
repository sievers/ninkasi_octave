function[crap]=ls_fixed(mystr)
[status,output]=system(['ls ' mystr ]);
if status==0,
  cc=char(10);
  crap=strsplit(strtrim(output),cc);
  if isempty(crap{end}),
    crap=crap(1:end-1);
  end
else
  crap={};
end

