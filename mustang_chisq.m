function[logl,mymod]=mustang_chisq(fitp)
global mynus
global mydat
mymod=fitp(1)+fitp(2)*mynus.^fitp(3);

logl=sum(mydat./mymod)+sum(log(mymod));

if ~isreal(logl)
  logl=1e8;
end


