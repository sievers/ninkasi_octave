function[matched]=match_tod_rowcol_pairs(tod,thresh)
if ~exist('thresh')
  thresh=0.1; %in arcminutes
end
[dx,dy]=get_detector_offsets_actpol(tod);
dx=repmat(dx,[1 numel(dx)]);
dy=repmat(dy,[1 numel(dy)]);
delt=(dx-dx').^2+(dy-dy').^2;
delt=sqrt(delt)*180/pi*60;
delt=delt+eye(length(dx))*max(max(delt));
[a,b]=min(delt);

matched=0*a;
%if distance is small enough to be considered a pair, find out which index it is
ismatched=a<thresh;
%now find which detector number we're at
matched(ismatched)=b(ismatched);

%finally, assign matched pairs.  shift by 1 to get to C indexing
set_detector_pairs_c(tod,matched-1);

%[rr,cc]