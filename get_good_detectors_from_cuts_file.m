function[rr,cc,cuts_org]=get_good_detectors_from_cuts_file(cutsname,maxnsamp)
if ~exist('maxnsamp')
  maxnsamp=1e9;
end

cuts_org=read_cuts_file_octave(cutsname);
cuts=cuts_org.det_cuts;
nrow=size(cuts,1);
ncol=size(cuts,2);
%[rr,cc]=meshgrid(1:size(cuts,1),1:size(cuts,2));
ok=false(nrow,ncol);
for j=1:nrow,
  for k=1:ncol
    vec=cuts{j,k};   
    if isempty(vec)
      ok(j,k)=true;
    else
      if max(isnan(vec))==0 %no nans present
        if (vec(1)~=0 ) | vec(2)<maxnsamp
          ok(j,k)=true;
        end
      end
    end      
  end
end
ngood=sum(sum(ok));
rr=zeros(ngood,1);
cc=zeros(ngood,1);
icur=0;
for j=1:nrow,
  for k=1:ncol,
    if ok(j,k),
      icur=icur+1;
      rr(icur)=j-1;
      cc(icur)=k-1;
    end
  end
end
disp(['have a total of ' num2str(sum(sum(ok))) ' OK detectors.']);


disp(['row/col is ' num2str([max(max(rr)) max(max(cc))])]);

