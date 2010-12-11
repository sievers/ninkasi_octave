function[dalt,dazcosalt]=convert_offsets_to_nomdelts(dx,dy,alt_ref)
%function[dalt,daz]=convert_offsets_to_nomdelts(dx,dy,alt_ref)
%dx,dy are in the (non-nomdelt) form, alt_ref is the median altitude of a TOD
bad_vals=abs(dx)>20;

dx(bad_vals)=0;
dy(bad_vals)=0;

dz=sqrt(1-dx.^2-dy.^2);
y2=dz*sin(alt_ref)+dy*cos(alt_ref);
z2=dz*cos(alt_ref)-dy*sin(alt_ref);
dalt=asin(y2)-alt_ref;
daz=atan2(dx,z2);
dazcosalt=daz*cos(alt_ref);


dalt(bad_vals)=-111;
dazcosalt(bad_vals)=-111;