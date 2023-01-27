function [LRs,rs_est,th_est] = test_Johansen_VECM_I2(y,p,alpha,const);
% I2_rank_det returns the reduced ranks, r and s, determinded by performing
% the LR rank test as well as the estimate obtained with this
% specification.
%
% SYNTAX:  [LRs,rs_est,th_est] = test_Johansen_VECM_I2(y,p,alpha,Joh_j);
%
% INPUT: y ... T x s real matrix of observations.
%        p ... integer, number of lagged differenced terms in VECM
%        alpha ... significance, \in {0.5, 0.2, 0.1, 0.05, 0.025, 0.01}.
%        Joh_j ... indicator: const = 1 means a constant is included.
%
% OUTPUT: LRs    ... (s+1) x (s+1) matrix of LRs. 
%         rs_est ... 3x1 integer of specified [r,s] values.
%         th_est ... theta strcture estimated using rs_est. 
%
% AUTHOR: Yuanyuan Li; adapted by dbauer, 26.8.2020. 
[T,s] = size(y);

if nargin<4
    const = 0;
end

if nargin <3
    alpha = 0.05;   % significance level
end;

% select the relevant critical values for inference.
cv = crit_val_sel(s,alpha,const);

r_hat = -1;
s_hat = -1;

for r = 0:s-1
    [th_est,S1,S2] = VECM_2SI2(y,p,[r,1],const);
    
    % [ r_hat, beta, alpha, beta_o, alpha_o, Omega, Q ]
    if r+1>length(S1)
        S1(r+1) = 0;
    end;
    
    lambda1 = S1.^2;        % (hat_lambda_1,...,hat_lambda_p)', descending.
    l_lambda1 = log(1-lambda1);
    Q1 =  - (T-p-2) * cumsum(l_lambda1, 'reverse');

    for d = 0:s-r
        if d+1>length(S2)
            S2(d+1) = 0;
        end;
        lambda2 = S2.^2;        % (hat_lambda_1,...,hat_lambda_p)', descending.
        l_lambda2 = log(1-lambda2);
        Q2 =  - (T-p-2) * cumsum(l_lambda2, 'reverse');
        Qrs = Q1(r+1) + Q2(d+1);
        LRs(r+1,d+1)= Qrs;
    end
end

for r = 0:s-1
    for d = 0:s-r
        Qrs = LRs(r+1,d+1);
        if Qrs <= cv(r+1,d+1)
            r_hat = r;
            s_hat = d;
            break;
        end        
    end
    
    if (r_hat>=0 || s_hat>=0)
        break;
    end;
    
end;

if (r_hat < 0 || s_hat < 0)
    r_hat = 0;
    s_hat = 0;
end;

rs_est = [r_hat,s_hat]; 
[th_est,Q1,Q2] = VECM_2SI2(y,p,rs_est,const);

