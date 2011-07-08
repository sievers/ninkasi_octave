function[kept_fracs]=cut_mostly_cut_detectors(tods,fac,nreg_thresh)
if length(tods)>1,
  kept_fracs={};
   for j=1:length(tods),
     kept_fracs(j)={cut_mostly_cut_detectors(tods(j))};
   end
   return;
end

[rows,cols]=get_tod_rowcol(tods);

nn=get_tod_ndata(tods);

if ~exist('fac')
  fac=0.5;
end
if ~exist('nreg_thresh')
  nreg_thresh=6;
end

kept_fracs=0*rows;
for j=1:length(rows),
    crud=print_tod_uncut_regions(tods,rows(j),cols(j));
    ll=crud(:,2)-crud(:,1);
    kept_fracs(j)=sum(ll)/nn;
    if (sum(ll))<fac*nn,
      cut_detector_c(tods,rows(j),cols(j));
      %      disp(['cutting ' [num2str(rows(j)) '  ' num2str(cols(j))] ' for being to short.']);
    end
    if length(ll)>nreg_thresh,
      %disp(['cutting ' [num2str(rows(j)) '  ' num2str(cols(j))] ' for being to chopped.']);
      cut_detector_c(tods,rows(j),cols(j));
    end
end

