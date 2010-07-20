function[lims]=get_tod_alldet_altaz_lims(tod)
if (1)
  lims=get_tod_alldet_altaz_lims_c(tod)';
  return
end

rows=get_tod_rowcol(tod);
ndet=length(rows);
lims=[1000 -1000 1000 -1000];
for j=0:ndet-1,
  altaz=get_detector_altaz_c(tod,j);
  altmin=min(altaz(:,1));
  altmax=max(altaz(:,1));
  azmin=min(altaz(:,2));
  azmax=max(altaz(:,2));
  if altmin<lims(1),
    lims(1)=altmin;
  end
  if altmax>lims(2)
    lims(2)=altmax;
  end
  if azmin<lims(3),
    lims(3)=azmin;
  end
  if azmax>lims(4),
    lims(4)=azmax;
  end
end
