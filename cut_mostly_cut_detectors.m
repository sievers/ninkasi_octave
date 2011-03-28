function[cut_fracs]=cut_mostly_cut_detectors(tods,fac)
if length(tods)>1,
  cut_fracs={};
   for j=1:length(tods),
     cut_fracs(j)={cut_mostly_cut_detectors(tods(j))};
   end
   return;
end

[rows,cols]=get_tod_rowcol(tods);

nn=get_tod_ndata(tods);

if ~exist('fac')
  fac=0.5;
end
nreg_thresh=6;
cut_facs=0*rows;
for j=1:length(rows),
    crud=print_tod_uncut_regions(tods,rows(j),cols(j));
    ll=crud(:,2)-crud(:,1);
    cut_fracs(j)=sum(ll)/nn;
    if (sum(ll))<fac*nn,
      cut_detector_c(tods,rows(j),cols(j));
      %      disp(['cutting ' [num2str(rows(j)) '  ' num2str(cols(j))] ' for being to short.']);
    end
    if length(ll)>nreg_thresh,
      %disp(['cutting ' [num2str(rows(j)) '  ' num2str(cols(j))] ' for being to chopped.']);
      cut_detector_c(tods,rows(j),cols(j));
    end
end

