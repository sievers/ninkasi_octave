function[nn]=get_tod_nuvec(tod)
dt=get_tod_dt(tod);
n=get_tod_ndata(tod);
nn=(0:n-1)';

nn=nn/(dt*n);
