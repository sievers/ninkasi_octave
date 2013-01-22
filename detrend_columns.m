function[mat]=detrend_columns(mat,npad)
if ~exist('npad')
  vec=(1:size(mat,1))'-1;vec=vec/(size(mat,1)-1);
  mat=mat-(repmat(mat(1,:),[size(mat,1) 1]) +  vec*(mat(end,:)-mat(1,:)));
else
  start_vals=mean(mat(1:npad,:),1); 
  stop_vals=mean(mat(end-npad+1:end,:),1);

  %extrapolate the line across the mean to the end
  slopes=(stop_vals-start_vals)/(length(mat)-npad);
  start_vals=start_vals-(npad-1)/2*slopes;
  stop_vals=stop_vals+(npad-1)/2*slopes;


  mat=mat-repmat(start_vals,[size(mat,1) 1]);

  stop_vals=stop_vals-start_vals;
  vec=(1:size(mat,1))'-1;vec=vec/(size(mat,1)-1);

  mat=mat-vec*stop_vals;
end

