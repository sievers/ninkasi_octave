function[map]=simple_allocate_ninkasi_skymap(pixsize,lims,varargin)
%simplified/more concise wrapper to allocate_ninkasi_skymap.  Args are (pixsize,lims,...) 
%where ... can include key/value pairs, such as 'pad' and the # of pixels by which to pad the map.
npad=get_keyval_default('pad',5,varargin{:});
pp=npad*pixsize;
map=allocate_ninkasi_skymap(pixsize,lims(1)-pp,lims(2)+pp,lims(3)-pp,lims(4)+pp);


