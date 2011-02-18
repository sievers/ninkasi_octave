function[lims,all_lims]=get_ra_wrapped_lims(all_lims,ra_wrap)
%shift TOD limits returned by read_all_tod_headers to be less than ra_wrap
ind =all_lims(:,1)>ra_wrap;
all_lims(ind,1)=all_lims(ind,1)-2*pi;
all_lims(ind,2)=all_lims(ind,2)-2*pi;

lims = [min(all_lims(:,1)) max(all_lims(:,2)) min(all_lims(:,3)) max(all_lims(:,4))];