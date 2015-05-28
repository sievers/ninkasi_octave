function[isrot]=are_dets_rotated(tod)

asdf=print_detector_pairs_c (tod)+1;ii=(1:length(asdf))';
%isrot=ii>asdf;
%isrot=~isrot;

isrot=(ii>asdf)&(asdf>0);