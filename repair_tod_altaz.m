function[alt,az,enc_flags]=repair_tod_altaz(tod)
if length(tod)>1,
  for j=1:length(tod),
    repair_tod_altaz(tod(j));
  end
  return
end

fname=get_tod_name(tod);
enc_flags=read_dirfile_channel(fname,'enc_flags','f');
if sum(enc_flags>0)==0  %nothing to see here, folks, move along.
  return;
end

alt=read_dirfile_channel(fname,'Enc_El_Deg','f');
az=read_dirfile_channel(fname,'Enc_Az_Deg','f');
alt=alt*pi/180;
az=(az+180)*pi/180;

inds=find(enc_flags~=0);

segs=inds2segs(inds')';

if size(segs,2)==1,
  assert(size(segs,1)==2);  %make sure it's really just a transpose issue
  segs=segs';
end


istart=segs(:,1)-1;
istop=segs(:,2)+1;


if istart(1)<1,
  istart(1)=istop(1);
end
if istop(end)>length(az)
  istop(end)=istart(end);
end

for j=1:size(segs,1),
  az(segs(j,1):segs(j,2))=myinterp1([segs(j,1)-1 segs(j,2)+1],[az(istart(j)) az(istop(j))],segs(j,1):segs(j,2));
  alt(segs(j,1):segs(j,2))=myinterp1([segs(j,1)-1 segs(j,2)+1],[alt(istart(j)) alt(istop(j))],segs(j,1):segs(j,2));
end

nn=get_tod_ndata(tod);
while length(az)>nn
  az=decimate_vector(az);
  alt=decimate_vector(alt);
end
assert(length(az)==nn);
push_tod_altaz(tod,alt,az);


return

function[yi]=myinterp1(x,y,xi)

dx=x(2)-x(1);dy=y(2)-y(1);yi=y(1)+dy*(xi-x(1))/dx;


