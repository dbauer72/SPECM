function test = eta_crit(s,j,T,M);
% eta_crit simulates the critical values for Johansen trace test. 
% 
% SYNTAX: test = eta_crit(s,j,T,M);
%
% INPUTS:   s ... integer; dimension of endogenous vars.
%           j ... integer: Johansens cases of deterministics: j in
%                (1,2,3,4,5).
%           T ... integer; sample size.
%           M ... integer; number of replications.
%
% OUTPUT: test ... M x s matrix of simulated values.
%
% REMARK: + test(:,k) contains the simulated values ordered in size for
%           dimension equal to k. 
%
% AUTHOR: dbauer, 14.1.2020. 


for m=1:M
    u = randn(T,s);
    y = cumsum(u);
    dy = y(2:end,:)-y(1:end-1,:);
    yl = y(1:end-1,:);
    pone = 0;
    switch j
        case 1
            yh=[1:T-1]'/T;
            yh=yh.^2;
            yl = [yh,yl];
            yl = detrend(yl,'linear');           
        case 2 % restricted linear trend 
            yh=[1:T-1]'/T;
            yl = [yh,yl];
            yl = detrend(yl,0);
            pone=1;
        case 3 % unrestricted constant
            yh=[1:T-1]'/T;
            yl = [yh,yl];
            yl = detrend(yl,0);
        case 4 % restricted constant.
            pone = 1;
            yh=ones(T-1,1);
            yl = [yh,yl];
        otherwise
            dy = y(2:end,:)-y(1:end-1,:);
            yl = y(1:end-1,:);
    end
    
    R00 = dy'*dy;
    R01 = dy'*yl;
    R11 = yl'*yl;
    
    for k=1:s % increase dimension
       inl = [1:(k+pone)];
       in = [1:k];
       test(m,k)= T*trace(inv(R11(inl,inl))*R01(in,inl)'*inv(R00(in,in))*R01(in,inl));
    end;
end;

for k=1:s
    test(:,k)=sort(test(:,k));
end;