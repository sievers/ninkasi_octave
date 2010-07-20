function[ra,dec]=evaluate_patchy_pointing_fit(alt,az,ctime,altvec,azvec,ra_mat,dec_mat,ra_clock,dec_clock)
ra=interp2(altvec',azvec,ra_mat',alt,az);
dec=interp2(altvec',azvec,dec_mat',alt,az);
if (1)
  ra=ra+polyval(ra_clock,ctime);
  dec=dec+polyval(dec_clock,ctime);
else
  ra=ra+ra_clock*ctime;
  dec=dec+dec_clock*ctime;
end
