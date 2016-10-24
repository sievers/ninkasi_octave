function[asdf]=harmonize_tod_detectors(tod1,tod2)
%cut all detectors that appear in only one of tod1/tod2
[rr1,cc1]=get_tod_rowcol(tod1);
[rr2,cc2]=get_tod_rowcol(tod2);
rc1=rr1+i*cc1;
rc2=rr2+i*cc2;
[ind,aa,bb]=setxor(rc1,rc2);
for j=1:length(aa),
  cut_detector_c(tod1,rr1(aa(j)),cc1(aa(j)));
end
for j=1:length(bb),
  cut_detector_c(tod2,rr2(bb(j)),cc2(bb(j)));
end
