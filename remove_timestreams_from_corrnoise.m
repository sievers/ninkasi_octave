function[new_corrnoise]=remove_timestreams_from_corrnoise(corrnoise,timestreams)
for j=size(corrnoise,2):-1:1,
  p=linfit(corrnoise(:,j),timestreams);
  new_corrnoise(:,j)=corrnoise(:,j)-timestreams*p;
end

