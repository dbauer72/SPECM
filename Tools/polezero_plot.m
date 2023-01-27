function ok = polezero_plot(ths,labs);
% polezero plots provides a plot of poles and zeros of a 
% number of transfer functions in state space form.
%
% SYNTAX: ok = polezero_plot(ths,labs);
%
% INPUT: ths ... vector of theta structures
%        labs ... vector of labels (used in the legend).
%
% OUTPUT: plot.
%
% AUTHOR: dbauer 25.3.2020

col = 'krbmgy';

k = length(ths);
if nargin< 2
    labs = cell(k,1);
    for j=1:k
        labs{j}= sprintf('th%d',j);
    end;
end;


fr = [0:1000]*2*pi/1000;
om = exp(sqrt(-1)*fr);
figure;
plot(real(om),imag(om));
hold on;

for j=1:k
    if isa(ths(1),'est_result')
        th = ths(j).theta;
    else
        th = poly2ss_th(ths(j));
    end
    
    A = th.A; 
    K = th.K;
    C = th.C;
    
    pol = eig(A);
    zero = eig(A-K*C);
    
    plot(real(pol), imag(pol),['x',col(rem(j,6)+1)]);
    plot(real(zero),imag(zero),['o',col(rem(j,6)+1)]);
end

labs_full = cell(0);
labs_full{1} = 'unit circle';
labs_full{2} = [labs{1}, ' Pole'];
labs_full{3} = [labs{1}, ' Zero'];
for j=2:k 
    labs_full{end+1}= [labs{j}, ' Pole'];
    labs_full{end+1}= [labs{j}, ' Zero'];
end;
    
legend(labs_full);