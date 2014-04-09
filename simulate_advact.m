function[asdf]=simulate_advact(tod,varargin)

knee=get_keyval_default('knee',1e-2,varargin{:});
powlaw=get_keyval_default('powlaw',1e-2,varargin{:});

dt=get_tod_dt(tod);
ndata=get_tod_ndata(tod);
nuvec=(0:ndata-1)';dnu=1/(dt*ndata);nuvec=nuvec*dnu;
dnu
