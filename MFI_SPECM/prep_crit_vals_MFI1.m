cd C:\dbauer\MATLAB\Projects\EICIP\SSECM\MFI_SSECM

% generate real valued critical values without deterministics
for n=1:11
    n
    [t_c]= simu_critval_MFI1(1000,100000,0,n,0);
    crit_val_real_wodet(n) = prctile(t_c,95);
end

% generate real valued critical values with deterministics
for n=1:11
    n
    [t_c]= simu_critval_MFI1(1000,100000,0,n,1);
    crit_val_real_det(n) = prctile(t_c,95);
end

% generate complex valued critical values without deterministics
for n=1:11
    n
    [t_c]= simu_critval_MFI1(1000,100000,0.25,n,0);
    crit_val_comp_wodet(n) = prctile(t_c,95);
end
% generate complex valued critical values with deterministics
for n=1:11
    n
    [t_c]= simu_critval_MFI1(1000,100000,0.25,n,1);
    crit_val_comp_det(n) = prctile(t_c,95);
end

save crit_MFI_ev_test crit_val_comp_det crit_val_comp_wodet crit_val_real_det crit_val_real_wodet