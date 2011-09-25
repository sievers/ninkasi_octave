function[taus]=assign_tod_time_constants(tods,taus)
if ~exist('taus')
  my_ar=get_tod_array(tods(1));assert(~isempty(my_ar));
  f3db=load(['/home/sievers/act/detectors/2008/' my_ar  '_f3db_090423.txt']);
  taus=1./f3db/2/pi;
  taus(f3db==0)=0;
  assert(sum(sum(isnan(taus)))==0);
  assert(sum(sum(isinf(taus)))==0);
else
  if ischar(taus)
    f3db=load(taus);
    taus=1./f3db/2/pi;
    taus(f3db==0)=0;
    assert(sum(sum(isnan(taus)))==0);
    assert(sum(sum(isinf(taus)))==0);
  else
    if iscell(taus)
      mdisp('assigning time constants from cell.');
      for j=1:length(tods),
        f3db=taus{guess_tod_season(tods(j))};
        mytau=1./f3db/2/pi;
        mytau(f3db==0)=0;
        assert(sum(sum(isnan(mytau)))==0);
        assert(sum(sum(isinf(mytau)))==0);
        assign_tod_time_constants_c(tods(j),mytau);
      end
      return;
    end
  end
end


for j=1:length(tods),
  assign_tod_time_constants_c(tods(j),taus);
end



 