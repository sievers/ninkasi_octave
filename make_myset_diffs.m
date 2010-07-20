function[maps,headers,names]=make_myset_diffs(fname,tag)
if ~exist('tag')
  tag='myset_';
end

pos=strfind(fname,tag);
if isempty(pos)
  error(['missing tag ' tag ' in file ' fname]);
end
ii=pos+length(tag);
myhead=fname(1:ii-1);
mytail=fname(ii+1:end);

imax=0;
while (1)
  ff=[myhead num2str(imax) mytail];
  fid=fopen(ff,'r');
  if (fid==-1)
    break;
  else
    fclose(fid);
    wname=get_weights_name(ff);
    fid=fopen(wname,'r');
    assert(fid~=-1);
    fclose(fid);
    imax=imax+1;
  end
end
disp(['have ' num2str(imax) ' sets.']);
assert(imax==4);

maps={};
headers={};
names={};
weights={};
for myset=1:imax,
  ff=[myhead num2str(myset-1) mytail];
  [map,header,nm]=fits_image_read(ff);
  maps(myset)={map};
  headers(myset)={header};
  names(myset)={nm};
  wname=get_weights_name(ff);
  ww=fits_image_read(wname);
  weights(myset)={ww};
end
whos

for ii=1:imax,
  for jj=ii+1:imax,
    map=maps{ii}-maps{jj};
    ww=1./(1./weights{ii}+1./weights{jj});
    ww(weights{ii}==0)=0;
    ww(weights{jj}==0)=0;
    outname=[myhead 'diff_' num2str(ii-1) '_' num2str(jj-1) mytail];
    write_fits(map,outname,names{ii},headers{ii});
    wname=get_weights_name(outname)
    write_fits(ww,wname,names{ii},headers{ii});
  end
end






function[wname]=get_weights_name(fname)
pos=strfind(fname,'.fits');
assert(pos>0);
pp=pos;
while fname(pp)~='_'
  pp=pp-1;
end
wname=[fname(1:pp) 'weights.fits'];
return