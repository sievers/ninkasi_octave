function[value]=write_corrnoise(corrnoise,outroot,tods)
if exist('tods')
  assert(length(corrnoise)==length(tods));
end
for j=1:length(corrnoise),
  fid=fopen([outroot '_' num2str(j) '.dat'],'w');
  if exist('tods')
    fwrite(fid,-1,'int');  %flag to say we have a more complete dump
    [rows,cols]=get_tod_rowcol(tods(j));
    fwrite(fid,length(rows),'int');
    fwrite(fid,rows,'int');
    fwrite(fid,cols,'int');
  end
    
  fwrite(fid,size(corrnoise(j).vecs),'int');

  for k=1:size(corrnoise(j).vecs,1),
    fwrite(fid,corrnoise(j).vecs(k,:),'double');
  end
  fwrite(fid,size(corrnoise(j).map),'int');
  for k=1:size(corrnoise(j).map,2),
    fwrite(fid,corrnoise(j).map(:,k),'double');
  end
  
  fclose(fid);
end
return

    