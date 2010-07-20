function[ra,dec]=get_pos_list_from_text_file(fname,varargin)
lines=read_text_file_comments(fname,varargin{:});
n=length(lines);
ra=zeros(n,1);
dec=zeros(n,1);
for j=1:n,
  [ra(j),dec(j)]=radecstring2deg(lines{j});
end


