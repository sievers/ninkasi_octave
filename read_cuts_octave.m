function[n_samp,samp_offset]=read_cuts_octave(tod,cutsname)
%C code was seg faulting on seemingly OK input, it was fragile to the presence of whitespace.  this should fix that...
ll=read_text_file(cutsname);
tags=strsplit(ll{1},' =',true)

if numel(tags)>=2
  if (strcmp(tags{1},'format'))&(strcmp(tags{2},'''TODCuts'''))
    [n_samp,samp_offset]=parse_cuts_octave_format(tod,ll);
    return
  end
end

[rr,cc]=get_tod_rowcol(tod);
max_rowcol=sscanf(ll{1},'%d');
assert(max(rr)<max_rowcol(1));
assert(max(cc)<max_rowcol(2));


found_cuts=false(size(rr));
nmax=get_tod_ndata(tod)

strip_tag=' ()rc:,';  %these characters need to die...
global_cuts=cellstr2mat(strsplit(ll{2},strip_tag,true));
if ~isempty(global_cuts)
  mdisp(['Have some global cuts in ' cutsname]);
  for j=1:2:length(global_cuts),
    for k=1:length(rr),
      cuts_extend_c(tod,global_cuts(j),global_cuts(j+1),rr(k),cc(k));
    end
  end
end

for j=3:length(ll)
  det_cuts=cellstr2mat(strsplit(ll{j},strip_tag,true));
  ii=find((det_cuts(1)==rr)&(det_cuts(2)==cc));
  have_det=~isempty(ii);
  found_cuts(ii)=true;
  myrow=det_cuts(1);
  mycol=det_cuts(2);
  mycuts=det_cuts(3:end);
  if (have_det) & numel(mycuts>0) %otherwise detector isn't cut at all



    %disp([myrow mycol])
    if (mycuts(1)==0)&&(mycuts(2)>nmax)
      %%disp('cutting detector entirely');
      cut_detector_c(tod,myrow,mycol);
    else
      for k=1:2:length(mycuts),
        cuts_extend_c(tod,mycuts(k),mycuts(k+1),myrow,mycol);
      end
    end
  end
end

if any(found_cuts==false)==false,
  %whos
  [rr,cc]=get_tod_rowcol(tod);
  for j=1:length(found_cuts),
    if found_cuts(j)==false,
      cut_detector_c(tod,rr(j),cc(j));
    end
  end
  %assert(any(found_cuts==false)==false)
end



