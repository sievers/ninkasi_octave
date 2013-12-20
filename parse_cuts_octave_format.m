function[n_samp,samp_offset]=parse_cuts_octave_format(tod,lines,varargin)
max_prime=get_keyval_default('max_prime',0,varargin{:});
j=1;
while (strcmp(strtrim(lines{j}),'END')==0)
  eval([lines{j} ';']);
  j=j+1;
end
if max_prime>0,
  facs=find_good_fft_lens(n_samp,max_prime);
  n_samp_good=max(facs);
  mdisp(['going from ' num2str(n_samp) ' to ' num2str(n_samp_good) ' samples.'])
  n_samp=n_samp_good;
end


%if we don't have a TOD, you can just get the offsets here without having to actually put in the cuts.
if isempty(tod)
  return;
end

jstart=j+1;
strip_tag=' rc:(),';

[rr,cc]=get_tod_rowcol(tod);
found_cuts=false(size(rr));
for j=jstart:length(lines),
  det_cuts=cellstr2mat(strsplit(lines{j},strip_tag,true));
  myind=find((det_cuts(2)==rr)&(det_cuts(3)==cc));
  if ~isempty(myind)
    found_cuts(myind)=true;
    myrow=det_cuts(2);
    mycol=det_cuts(3);
    if length(det_cuts)>3 %if nothing, it isn't cut
      mycuts=det_cuts(4:end);
      if (mycuts(1)==0)&&(mycuts(2)>=n_samp)
        cut_detector_c(tod,myrow,mycol);
      else
        for k=1:2:length(mycuts),
          cuts_extend_c(tod,mycuts(k),mycuts(k+1),myrow,mycol);
        end
      end
    end
  end

end
rr_bad=rr(found_cuts==false);
cc_bad=cc(found_cuts==false);
if (~isempty(rr_bad))
  mdisp(['cutting ' num2str(length(rr_bad)) ' detectors for missing cuts.'])
  for j=1:length(rr_bad),
    cut_detector_c(tod,rr_bad(j),cc_bad(j));
  end
end

