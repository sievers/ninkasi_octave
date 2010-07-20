function[value]=plot_detectors(data,rows,cols,varargin)
if length(varargin)==1
    ind=varargin{1};
else
    myrow=varargin{1};
    mycol=varargin{2};
    ii=1:size(data,2);
    myj=0;
    for j=1:length(myrow),
        jj=(myrow(j)==rows)&(mycol(j)==cols);
        if (sum(jj)==1)
            myj=myj+1;
            ind(myj)=ii(jj);
        end
    end
end


my_legend={};
for j=1:length(ind),
    rr=rows(ind(j));
    cc=cols(ind(j));
    ff=['r:' num2str(rr) ' c:' num2str(cc)];
    my_legend(end+1)={ff};
end


crap=data(:,ind);
dx=mean(mean(abs(crap)));
cm=median(data,2);
crap(:,end+1)=cm;
for j=1:size(crap,2),
    crap(:,j)=crap(:,j)+4*j*dx;
end

my_legend(end+1)={'common mode'};
plot(crap)
legend(my_legend{:},'location','best');
