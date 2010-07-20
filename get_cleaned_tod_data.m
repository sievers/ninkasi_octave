function[x]=get_cleaned_tod_data(tod)
read_tod_data(tod);
gapfill_data_c(tods);
detrend_data(tods);
window_data(tods);
