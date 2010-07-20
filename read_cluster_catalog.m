function[mycat]=read_cluster_catalog(fname,fdir)
if ~exist('fdir')
  fdir='/home/sievers/act/catalogs/';
end
if ~exist('fname')
  fname='final_catalog_for_zooms.txt';
end

if sum(fname=='/')==0
  fname=[fdir '/' fname];
end
ll=read_text_file(fname);
if isempty(ll)
  error(['File ' fname ' was empty in read_cluster_catalog.']);
end
ncur=0;
n=length(ll);
mycat(n).name='';
mycat(n).ra=0;
mycat(n).dec=0;
mycat(n).extras='';
for j=1:length(ll)
  a=strtrim(ll{j});
  if sum(a=='#')>0
    a=a(1:min(find(a=='#'))-1);
  end
  if ~isempty(a)
    ncur=ncur+1;
    [nm,a]=strtok(a);
    [ra,a]=strtok(a);
    [dec,a]=strtok(a);
    [snr,extras]=strtok(a);
    mycat(ncur).name=nm;
    mycat(ncur).ra=str2num(ra);
    mycat(ncur).dec=str2num(dec);
    mycat(ncur).extras=strtrim(extras);
  end
end

mycat=mycat(1:ncur);
  
