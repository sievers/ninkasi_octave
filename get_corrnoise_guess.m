function[corrnoise]=get_corrnoise_guess(data)

corrmat=data'*data;

[v,e]=eig(corrmat);
e=diag(e);
thresh=10;  %treat any modes with eigenvalue more than this times the median as special
v_use=v(:,e>thresh*median(abs(e)));
s_guess=data*v_use;

corrnoise.map=s_guess;
corrnoise.vecs=v_use';
