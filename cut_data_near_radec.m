function[value]=cut_data_near_radec(tod,ra,dec,radius)
if length(tod)>1
  for j=1:length(tod),
    cut_data_near_radec(tod(j),ra,dec,radius);
  end
  return
end

[rows,cols]=get_tod_rowcol(tod);
if ra>pi,
  ra2=ra-2*pi;
else
  ra2=ra+2*pi;
end

for j=1:length(rows),
  radec=get_detector_radec(tod,rows(j),cols(j));
  radec(radec(:,1)<0,1)=radec(radec(:,1)<0,1)+2*pi;
  

  %dist=(ra-radec(:,1)).^2*cos(dec)^2+(radec(:,2)-dec).^2;
  dra1=(ra-radec(:,1)).^2;
  dra2=(ra2-radec(:,1)).^2;
  dist=min(dra1,dra2)*cos(dec)^2+(radec(:,2)-dec).^2;
  crud=diff(dist<radius^2);
  istart=find(crud==1);
  istop=find(crud==-1);

  if ~isempty(istart)
    if istop(1)<istart(1)
      istart=[1 istart];
    end
    if istop(end)<istart(end)
      istop=[istop length(dist)];
    end
    assert(length(istart)==length(istop));
    for k=1:length(istart),
      cuts_extend_c(tod,istart(k),istop(k),rows(j),cols(j));
    end
  end
end


