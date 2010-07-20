function[value]=get_rsync_stuff(tod_names)
for j=1:length(tod_names),

    tt=tod_names{j};
    tt2=tt;tt2(tt2=='.')=' ';tt2=str2num(strtok(tt2));
    a=date_from_ctime(tt2,'hh');
    a=str2num(a);
    dd=date_from_ctime(tt2,'yyyymmdd');
    if a>12
       dd=str2num(dd);
       dd=dd+1;
       dd=sprintf('%d',dd);
     end
     disp([dd '/' tod_names{j}]);

end;
