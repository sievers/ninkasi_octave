function[vecs]=find_bad_modes_block(datft,nu,freqs,eig_thresh,skip_mean)

if (~skip_mean)
  vecs=ones(size(datft,2),1);  %explicitly project the mean
  for j=1:size(vecs,2),
    vecs(:,j)=vecs(:,j)/sqrt(sum(vecs(:,j).^2));
  end
else
  vecs=[];
end


for j=1:length(freqs)-1,
  ind=(nu>freqs(j))&(nu<=freqs(j+1));
  ind_lims=[min(find(ind)) max(find(ind))];
  %mdisp(num2str(ind_lims));
  crud=datft(ind,:);

  mat=real(crud'*crud);
  mat=mat+mat';
  if (~isempty(vecs))
    mat=project_vecs_from_mat(mat,vecs);
  end
  mat=0.5*(mat+mat');

  if sum(sum(~isfinite(mat)))>0
    %crap=get_tod_name(tod);
    %while iscell(crap)
    %  crap=crap{1};
    %end
    %disp('we have a problem in find_bad_mode opts: ')
    %disp(crap)
    vecs=[];
    return
  end
  [vv,ee]=eig(mat);
  ee=diag(ee);
  %put in a check so that if our frequency bin is too narrow we don't puke
  nind=ind_lims(2)-ind_lims(1)+1;  
  if nind<length(ee)
    mdisp(['not enough modes in band 1, using ' num2str(nind) ' modes in find_bad_modes_opts.m']);
    ee_med=median(ee(end-nind+1:end));
  else
    ee_med=median(ee);
  end

  if eig_thresh(j)>0,
    val=eig_thresh(j)*ee_med;
    ind=ee>val;
  else
    ind=false(size(ee));
    ind(end+eig_thresh+1:end)=true;
  end
  vecs=[vecs vv(:,ind)];
    
  mdisp(sprintf('cutting %d out of %d modes on band %d',sum(ind),length(ind),j));
  clear crud
end

