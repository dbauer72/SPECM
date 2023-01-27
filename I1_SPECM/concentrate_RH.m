function [R0t,R1t] = concentrate_RH(y,Abar,B,rest,Joh_j); 
% Concentration step for Ribarits_Hanzon idea for estimation of RH_SPECM form
%
% SYNTAX: [R0t,R1t] = concentrate_RH(y,Abar,B,rest,Joh_j); 
%
% INPUT: y    ... Txs data matrix,
%        Abar ... nxn matrix A-BC.
%        B    ... nxs matrix
%        rest ... 'y'|'n' restricted version.
%        Joh_j ... integer; indicating the version of Johansen's
%        deterministics. 
%
% OUTPUT: R0t, R1t ... Txs matrices of concentrated processes Z0t and Z1t. 
%
% AUTHOR: dbauer, 8.2.2022. 

[T,s] = size(y);
n = size(Abar,1);
% steps: calculate tilde x_t: filtered version of Delta y.

Z0t = diff(y);Z0t = detrend(Z0t,0);
Z1t = y(1:end-1,:);

Z2t = zeros(T,n);

tB = B';
tA = Abar';
for t=2:T
    Z2t(t,:)=Z2t(t-1,:)*tA + Z0t(t-1,:)*tB;
end;
%Z2t = ltitr(Abar,B,Z0t,zeros(n,1));

Z2t = Z2t(2:end-1,:)*inv(eye(n)-tA)*tA;
Z0t = Z0t(2:end,:);
Z1t = Z1t(2:end,:);

plus1 = 0;
switch Joh_j
    case 1
        Z0t = detrend(Z0t,1);
        Z1t = detrend(Z1t,1);
        Z2t = detrend(Z2t,1);
    case 2 % restricted linear trend
        Z0t = detrend(Z0t,0);
        Z1t = [[1:T-2]',Z1t];
        Z1t = detrend(Z1t,0);
        Z2t = detrend(Z2t,0);
        plus1=1;
    case 3 % unrestricted constant
        Z0t = detrend(Z0t,0);
        Z1t = detrend(Z1t,0);
        Z2t = detrend(Z2t,0);
    case 4
        Z1t = [ones(T-2,1),Z1t];
        plus1=1;
end

switch rest
    case 'y'
        % set up restriction
        R = inv(eye(n)-Abar)*B;
        
        % add restrictions, if not enough
        if n<s % less states than output dimensions
            [LL,alphahat,betahat,Chat] = RH_vecm(y,Abar,B,r,'n');
            beta2 = betahat(:,1:(s-n));
            R = [R;beta2'];
            Z2t = [Z2t,Z1t*beta2];
            r = r-s+n;
        end;
        % regressor matrix X'X:
        tildeM = Z2t'*Z2t/T;
        itM = inv(tildeM);
        M02 = Z0t'*Z2t/T;
        M12 = Z1t'*Z2t/T;
        
        C0hat = M02*itM;
        C1hat = M12*itM;
        
        H11 = itM*(eye(size(R,1))-R*inv(R'*itM*R)*R'*itM);
        H21 = inv(R'*itM*R)*R'*itM;
        
        C0bar = M02*H11 - H21;
        C1bar = M12*H11 + [zeros(plus1,size(H21,2));H21];
        
        % generate R_{0t}
        R0t = Z0t - Z2t*C0bar';
        R1t = Z1t - Z2t*C1bar';
    otherwise % here comes the unrestricted estimation.
        C0 = (Z2t\Z0t);
        C1 = (Z2t\Z1t);
        
        % generate R_{0t}
        R0t = Z0t - Z2t*C0;
        R1t = Z1t - Z2t*C1;
end

