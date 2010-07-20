function[value]=highpass_tod(tods,t1,t2);
if ~exist('t1','var')
  t1=0.5;
end
if ~exist('t2','var')
  t2=2*t1;
end

for j=1:length(tods),
  highpass_tod_c(tods(j),t1,t2);
end
