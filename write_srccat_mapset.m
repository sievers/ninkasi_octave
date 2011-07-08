function[value]=write_srccat_mapset(srccat,outfile)
if iscell(srccat),
  for j=1:numel(srccat),
    write_srccat_mapset(srccat{j},[outfile '._cat_' num2str(j)]);
  end
  return
end

fid=fopen(outfile,'w');
ra=srccat.ra;
dec=srccat.dec;
amps=srccat.amps;
nsrc=numel(ra);
for j=1:nsrc,
  if (1)
    fprintf(fid,'%12.6f %12.6f %15.6e\n',ra(j)*180/pi,dec(j)*180/pi,amps(j));
  else
    
    rr=ra(j)*180/pi/15;
    rah=floor(rr);
    rr=(rr-rah)*60;
    ram=floor(rr);
    ras=60*(rr-ram);
    
    dd=dec(j)*180/pi;
    isign=dd>0;
    dd=abs(dd);
    ddeg=floor(dd);
    dd=(dd-ddeg)*60;
    dm=floor(dd);
    ds=60*(dd-dm);
    
    if isign>0
      signchar='+';
    else
      signchar='-';
    end
    fprintf(fid,'%2d %2d %6.3f %c%0d %2d %5.2f %14.5g\n',rah,ram,ras,signchar,ddeg,dm,ds,amps(j));
  end
end
fclose(fid);
