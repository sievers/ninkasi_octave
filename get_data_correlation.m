function[corrs]=get_data_correlation(mytod)
corrs=get_data_correlation_c(mytod);
corrs=corrs+corrs';corrs=corrs-diag(0.5*diag(corrs));


