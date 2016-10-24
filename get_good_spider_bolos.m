function[asdf]=get_good_spider_bolos(fname,focal_plane)
if ~exist('fname')
  fname='/home/cynthia/analysis/spider_good_bolos.txt';
end

asdf=parse_spider_offsets(fname);

if exist('focal_plane')
  asdf=asdf(asdf(:,1)==focal_plane,:);
end

