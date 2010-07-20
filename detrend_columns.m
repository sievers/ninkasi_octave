function[mat]=detrend_columns(mat)
vec=(1:size(mat,1))'-1;vec=vec/(size(mat,1)-1);
mat=mat-(repmat(mat(1,:),[size(mat,1) 1]) +  vec*(mat(end,:)-mat(1,:)));


