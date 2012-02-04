function[myrow,mycol]=find_nearest_live_detector(tod,row,col)
[rows,cols]=get_tod_rowcol(tod);
dist=sqrt( (rows-row).^2+(cols-col).^2);
[a,b]=min(dist);
myrow=rows(b);
mycol=cols(b);
