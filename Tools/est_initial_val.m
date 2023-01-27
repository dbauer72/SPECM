function [res,x0_est] = est_initial_val(res,Abar,tilK,tilC);
% est_initial_val estimates the initial values for a state space system 
% in the PEM approach. 
%
% SYNTAX: [res,x0_est] = est_initial_val(res,Abar,tilK,tilC);
%
% INPUT: res ... Txs real matrix of residuals;
%        (Aber,tilK,tilC) state space system of the inverse transfer
%        function. 
%
% OUTPUT: res ... Txs real matrix of updated residuals.
%         x0_est ... nx1 estimated initial values.
%
% AUTHOR: dbauer, 8.10.2021.

[Ts,s]=size(res);
n = size(Abar,1);

% estimate initial state to get rid of initial effects
x0_all = eye(n);

res_v = res(:); % vectorization of residuals.
X = zeros(length(res_v),n);
for j=1:n
    xinj = ltitr(Abar,tilK,res(1:Ts,:)*0,x0_all(:,j))*tilC';
    X(:,j) = xinj(:);
end

x0_est = X\res_v;
tres_v = res_v - X*x0_est;

res = reshape(tres_v,Ts,s);

