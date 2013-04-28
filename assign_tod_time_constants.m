function[taus]=assign_tod_time_constants(tods,taus)

myid=mpi_comm_rank;

if ~exist('taus')
  my_ar=get_tod_array(tods(1));assert(~isempty(my_ar));
%f3db=load(['/home/sievers/act/detectors/2008/' my_ar  '_f3db_090423.txt']);
  f3db=load(['/home/r/rbond/sievers/act/detectors/2008/' my_ar  '_f3db_090423.txt']);
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
      %disp(num2str([mpi_comm_rank size(taus{4})]))
      mpi_barrier;
      
      mdisp('assigning time constants from cell.');
      for j=1:length(tods),
        %disp(['process ' num2str(myid) ' starting on ' num2str(j)]);
        f3db=taus{guess_tod_season(tods(j))};
        mytau=1./f3db/2/pi;
        mytau(f3db==0)=0;
        assert(sum(sum(isnan(mytau)))==0);
        assert(sum(sum(isinf(mytau)))==0);
        %mpi_barrier;mdisp(['j is ' num2str(j)]);
        assign_tod_time_constants_c(tods(j),mytau);
        %disp(['process ' num2str(myid) ' finished on ' num2str(j)]);
        %mpi_barrier;mdisp('passed.');
      end
      return;
    end
  end
end



for j=1:length(tods),
  assign_tod_time_constants_c(tods(j),taus);
end



 
