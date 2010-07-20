function[x,a,b]=simple_solve_mat_test(a,b)
%solve  ax=b

n=length(b);
for j=1:n-1,
  for k=j+1:n,
    fac=a(k,j)/a(j,j);
    a(k,:)=a(k,:)-a(j,:)*fac;
    b(k)=b(k)-b(j)*fac;
  end
end
x=0*b;
x(n)=b(n)/a(n,n);
for j=n-1:-1:1,
  x(j)=(b(j)-a(j,j+1:end)*x(j+1:end))/a(j,j);
end




    