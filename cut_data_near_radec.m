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
    if isempty(istop)
      istop=length(dist);
    end
    %JLS - fixing corner case where istop(1)==1 - 26May 2011
    if istop(1)<istart(1);
      if numel(istop)>1
        istop=istop(2:end)
      else
        istop=length(dist); %import empty case from above
      end      
    end
    %End mod
    try
      if istop(1)<istart(1)
        istart=[1 istart];
      end
    catch
      whos
      istart
      istop
      disp(['screwed up again. on ' get_tod_name(tod)]);
      error(['screwed up again. on ' get_tod_name(tod)]);
    end
    if istop(end)<istart(end)
      try
        %istop=[istop length(dist)];
        istop(end+1)=length(dist);
      catch
        whos
        istart
        istop
        error(['tod ' get_tod_name(tod) ' is screwed up.']);
        
        return
      end
    end
    assert(length(istart)==length(istop));
    for k=1:length(istart),
      cuts_extend_c(tod,istart(k),istop(k),rows(j),cols(j));
    end
  end
end


