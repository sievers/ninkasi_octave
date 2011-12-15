function[names,dists]=convert_namedist(lines,ra,dec)
names={};
dists=zeros(numel(lines),1);
for j=length(lines):-1:1,
  [nm,crud]=strtok(lines{j});
  mypos=str2num(crud);
  names(j)={nm};
  dra=(mypos(1:2)-ra)*cos(dec);
  ddec=(mypos(3:4)-dec);
  mat=repmat(dra.^2,[2 1])+repmat(ddec'.^2,[1 2]);
  dists(j)=min(min(mat));
end

dists=sqrt(dists);
