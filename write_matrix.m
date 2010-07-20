function[value]=write_matrix(mat,fname)

fid=fopen(fname,'w');
fwrite(fid,size(mat),'int');
fwrite(fid,mat,'double');
fclose(fid);

