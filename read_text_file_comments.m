function[lines]=read_text_file_comments(fname,comment_str)
if ~exist('comment_str')
  comment_str='#';
end
lines=read_text_file(fname);
for j=1:length(comment_str),
  for k=1:length(lines),
    crud=lines{k};
    fwee=find(crud==comment_str(j));
    if ~isempty(fwee)
      crud=crud(1:min(fwee)-1);
      lines(k)={crud};
    end
  end
end
ind=true(size(lines));
for j=1:length(lines),
  if isempty(strtrim(lines{j}))  %used to be just isempty - now blank lines w/spaces should get cut
    ind(j)=false;
  end
end
lines=lines(ind);

