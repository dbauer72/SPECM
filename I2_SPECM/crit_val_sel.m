function [ crit_val ] = crit_val_sel( p,alpha, const )
% For given a significance level, alpha, 'crit_val_sel' returns a matrix of
% critical valutes of the trace statistic for I2, Q_r+s, obtained from the
% 2-step RRR procedure.
% 
% SYNTAX: crit_val = crit_val_sel(p,alpha);
%
% INPUT : p ... integer, dimension of the data, [1,12].
%         alpha ...significance level, {0.5, 0.2, 0.1, 0.05, 0.025, 0.01}.
%
% OUTPUT: crit_val ... p x (p+1) matrix; crit_val(i,j): critical value of
%                (r=i-1, s1=j-1, s2=p-r-s1).
%
% EXTERNAL FUNCTIONS: 'cv_table2.mat'/I2_crit_val_calc(to be done).

if p>12 || p<=0
    error('For p=1,...,12.');
end;

crit_val = zeros(p,p+1);

sig_alpha = [0.5, 0.2, 0.1, 0.05, 0.025, 0.01];
if const == 0
    cv_file = matfile('cv_table2.mat');
    C = cv_file.cv_table2;
else 
    cv_file = matfile('cv_table2_c.mat');
    C = cv_file.cv_table2;
end

k = find(sig_alpha == alpha);
if isempty(k)
    error('alpha should be one of the following: 0.5, 0.2, 0.1, 0.05, 0.025, 0.01.');
end;

n = (3+p)*p/2;    % n = 2+3+...+(p+1).
cv_sel = C(end-n+1:end,k);

pp=1;
for i = 1:p
    crit_val(i,1:p+2-i) = cv_sel(pp:pp+p+1-i)';
    pp = pp+p+2-i;
end;

end

