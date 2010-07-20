function[pp,m1,m2]=decompose_2mats(a1,a2)
%calculate pp,m1,m2 such that a1=pp*m1*pp', a2=pp*m2*pp', m1, m2
%are diagonal, and pp has columns of length one.

[v1,ee]=eig(a1);

if min(min(ee))<0
  error('Please send in positive definite matrices.');
end


s=sqrt(ee);


bb=inv(s)*v1'*a2*v1*inv(s);
bb=(bb+bb')/2;
[u,m]=eig(bb);
u=u';
pp=v1*s*u';
pplen=(sum(pp.^2,1));
m1=pplen';
m2=diag(m).*(pplen');


pp=pp*diag((1./sqrt(pplen)));




