addpath /home/sievers/matlab
more off
tic;tod_names=read_text_file('/home/sievers/tods_rising_1.txt');toc
tod_names=[tod_names read_text_file('/home/sievers/tods_setting_1.txt')];toc
tod_names=[tod_names read_text_file('/home/sievers/tods_rising_0.txt')];toc
tod_names=[tod_names read_text_file('/home/sievers/tods_setting_0.txt')];toc

%tod_names=tod_names(1:100);

for j=length(tod_names):-1:1,
	c1{j}=guess_cuts_name(tod_names{j},'/home/sievers/act/cuts/season2_cuts/');
        c2{j}=guess_cuts_name(tod_names{j},'/home/sievers/act/cuts/season2_cuts3/');
	c3{j}=guess_cuts_name(tod_names{j},'/home/sievers/act/cuts/season2_cuts_cuts3_merged/');
end
merge_cuts_c(c1,c2,c3);
