function[ra,dec]=read_radec_from_txtfile(fname)
lines=read_text_file_comments(fname);
ra=zeros(length(lines),1);
dec=ra;
for j=1:length(lines),
  [ra(j),dec(j)]=parse_radec_text_line(lines{j});
end
