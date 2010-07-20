function[value,projected]=project_vecs_from_mat(mat,vecs)
if (1)
  vec=mat*vecs;
  a=(vecs'*vec);
  projected=vec*inv(a)*vec';
  value=mat-projected;
else
  [u,s,v]=svd(vecs);

  uu=u(:,size(vecs,2)+1:end);

  mm=uu'*mat*uu;
  value=uu*mm*uu';

  if (nargout>1)
    uu=u(:,1:size(vecs,2));
    mm=uu'*mat*uu;
    projected=uu*mm*uu';
  end
end
