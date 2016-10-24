function[vals]=read_cfg(fname)
ll=read_text_file_comments(fname);
ii=zeros(size(ll));
for j=1:length(ll),
  ii(j)=ll{j}(1)=='[';
end
ii=find(ii);
nb=length(ii);
ii(end+1)=length(ll)+1;

vals=cell(nb,1);
for j=1:nb,
  %disp(j)
  block=ll(ii(j):ii(j+1)-1);
  ln=strtrim(block{1});
  %disp(ln)
  clear mystruct
  mystruct.name=ln(2:end-1);
  for jj=2:length(block),
    aa=min(find(block{jj}=='='));
    tag=strtrim(block{jj}(1:aa-1));
    val=strtrim(block{jj}(aa+1:end));
    %if there are quotes in the field, fix them
    val=escape_quotes(val);
    to_eval=['mystruct.' tag ' = ''' val ''';'];
    eval(to_eval);
  end
  vals(j)=mystruct;

end


function[asdf]=escape_quotes(val)
%replace all single quotes with two single quotes so that the eval above handles them correctly
ii=find(val=='''');
if isempty(ii)
  asdf=val;
  return;
end

asdf=[val(1:ii(1)) ''''];
for j=2:length(ii),
  asdf=[asdf val(ii(j-1)+1:ii(j)) ''''];
end
if length(val)>ii(end),
  asdf=[asdf val(ii(end)+1:end)];
end
