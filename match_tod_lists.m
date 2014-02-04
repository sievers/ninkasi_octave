function[tods,vec]=match_tod_lists(targs,tod_list)
%find the tods in the second list that match the first.  Useful for, e.g., getting full path names for a list of tods by ctime
[ct,ct2]=get_tod_ctimes_from_names(tod_list);
[tt,tt2]=get_tod_ctimes_from_names(targs);
vec=0*tt;
for j=1:length(tt),
  ii=find(tt(j)==ct);
  if numel(ii)==1,
    vec(j)=ii;
  else
    ii2=find(tt2(j)==ct2);
    myind=intersect(ii,ii2);
    if numel(myind)==1
      vec(j)=myind;
    else
      vec(j)=nan;
    end
  end
end
if sum(isnan(vec)>0)
  mdisp(['Warning - was missing ' num2str(sum(isnan(vec))) ' matches in match_tod_lists.']);
end

tods=tod_list(vec(isfinite(vec)));

