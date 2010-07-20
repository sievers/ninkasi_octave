function[map2,names,vals,map]=read_map_highpassed(fname,varargin)
thresh=get_keyval_default('thresh',0.1,varargin{:});
width=get_keyval_default('width',0.2,varargin{:});

[map,names,vals]=fits_image_read(fname);

if thresh>0
  weight_name=guess_weight_name(fname);
  if (weight_name)    
    wt=fits_image_read(weight_name);
    wt=wt/(thresh*max(max(wt)));
    wt(wt>1)=1;
    map=map.*wt;
  end
end



mapft=fft2(map);



vec=(1:size(map,1));vec=vec-ceil((size(map,1)+1)/2);x=ifftshift(vec)';
vec=(1:size(map,2));vec=vec-ceil((size(map,2)+1)/2);y=ifftshift(vec);
x=x/max(x);
y=y/max(y);
mat=repmat(x.^2,[1 size(map,2)]);
mat=mat+repmat(y.^2,[size(map,1) 1]);
mapft=mapft.*(1-exp(-0.5*mat/width^2));
clear mat;
map2=real(ifft2(mapft));
whos






return


function[nm2]=guess_weight_name(nm)

ind=max(find(nm=='_'));
nm2=[nm(1:ind) 'weights.fits'];
fid=fopen(nm2,'r');
if fid<0
  nm2='';
else
  fclose(fid);
end
