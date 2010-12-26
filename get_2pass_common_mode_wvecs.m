function[cm,cm_amps,vec_amps,data_clean,cm_org]=get_2pass_common_mode_wvecs(data,vecs)
cm=median(data')';

vv=[vecs cm];
atd=vv'*data;
ata=vv'*vv;
fitp=inv(ata)*atd;
data_scale=data;for j=1:size(data,2),data_scale(:,j)=data_scale(:,j)/fitp(end,j);end;

%pull out the solved-for vectors from the rescaled_data;
atd=vecs'*data_scale;
ata=vecs'*vecs;
fitp=inv(ata)*atd;
data_scale=data_scale-vecs*fitp;


cm_org=cm;
cm=median(data_scale')';

%pull out the solved-for vecs from the common mode
atd=vecs'*cm;
ata=vecs'*vecs;
fitp=inv(ata)*atd;
cm=cm-vecs*fitp;



vv=[vecs cm];
atd=vv'*data;
ata=vv'*vv;
fitp=inv(ata)*atd;


data_clean=data-vv*fitp;
cm_amps=fitp(end,:);
vec_amps=fitp(1:end-1,:);

