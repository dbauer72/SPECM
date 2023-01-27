function [LL,alphahat,betahat,Chat] = RH_vecm_H0(y,Abar,B,r,rest,Joh_j,restrict);
% ribarits_hanzon idea for estimation of RH_VECM form
%
% SYNTAX: [LL,alphahat,betahat,Chat] = RH_vecm_H0(y,Abar,B,r,rest,Joh_j,restrict);
%
% INPUT: y    ... Txs data matrix,
%        Abar ... nxn matrix A-BC.
%        B    ... nxs matrix
%        r    ... rank of alpha beta' 
%        rest ... 'y'|'n' restricted version.
%        Joh_j ... integer; indicating the version of Johansen's
%        deterministics. 
%        restrict ... takes restrictions into account. restrict.H carries
%                  information.
%
% OUTPUT: LL   ... minimized log likelihood.
%         alphahat, betahat ... sxr matrices.
%         Chat ... sxn optimal C.
%
%
% COMMENT: + The restriction C*(I-Abar)^(-1)*B = I+alpha beta' is imposed by
% rest = 'y'.  
%          + implements estimation for betahat = H varphi, H specified in
%          restrict: restrict.H in R^(s x t), restrict.t = t.
%
% AUTHOR: dbauer, 21.11.2017. 

% prepare regressions 
if nargin< 6
    Joh_j = 5; % no deterministic terms. 
end;

if nargin <7
    H= [];
else
    if isfield(restrict,'H')
        H= restrict.H;
        tH = restrict.t;
    else
        H = [];
    end;
end;

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

Z2t = Z2t(2:end-1,:)*inv(eye(n)-tA)*tA;
Z0t = Z0t(2:end,:);
Z1t = Z1t(2:end,:);

switch Joh_j
    case 1
        Z0t = detrend(Z0t,'linear');
        Z1t = detrend(Z1t,'linear');
        Z2t = detrend(Z2t,'linear');
    case 2 % restricted linear trend
        Z0t = detrend(Z0t,0);
        Z1t = detrend(Z1t,0);
        Z2t = detrend(Z2t,0);
    case 3 % unrestricted constant
        Z0t = detrend(Z0t,0);
        Z1t = detrend(Z1t,0);
        Z2t = detrend(Z2t,0);
    otherwise
end

% here comes the separation corresponding to restriction or not.

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
        end;
        
        
        % differentiate, if columns are left:
        if r >= tH % no additional columns left
            if (r<tH)
                disp('RH_vecm_H0: Specification is wrong! r<t.');
                return;
            else % no r==t, restriction fully defines cointegrating space -> OLS regression.
                Z1th = Z1t*H(:,1:tH);
                % adjust for deterministics.
                plus1 = 0;
                switch Joh_j
                    case 2
                        yh=[1:(T-2)]'/(T-2);
                        Z1th = [yh,Z1th];
                        plus1 = 1;
                    case 4
                        yh=ones(T-2,1);
                        Z1th = [yh,Z1th];
                        plus1 = 1;
                end;
                
                % calculate regression
                betahat = H(:,1:tH);
                alph = ([Z1th,Z2t]\Z0t);
                alphahat = alph(plus1+[1:tH],:)';
                res = Z0t - [Z1th,Z2t]*alph;
                LL =  (T-2)*log(det(res'*res/(T-2)))+(T-2)*s;
                Chat = -alph(plus1+tH+1:end,:)';
            end
            
            return
        end
        
        % now t>r!
        Z1th = Z1t*H(:,1:tH);
        
        % regressor matrix X'X:
        tildeM = Z2t'*Z2t/T;
        itM = inv(tildeM);
        M02 = Z0t'*Z2t/T;
        M12 = Z1th'*Z2t/T;
        
        C0hat = M02*itM;
        C1hat = M12*itM;
        
        H11 = itM*(eye(n)-R*inv(R'*itM*R)*R'*itM);
        H21 = inv(R'*itM*R)*R'*itM;
        
        C0bar = M02*H11 - H21;
        C1bar = M12*H11 + H(:,1:tH)'*H21;
        %        C0bar = C0hat - (C0hat*R+eye(s))*inv(R'*itM*R)*R'*itM;
        %        C1bar = C1hat - (C1hat*R-eye(s))*inv(R'*itM*R)*R'*itM;
        
        % generate R_{0t}
        R0t = Z0t - Z2t*C0bar';
        
        % generate R_{1t}
        R1t = Z1th - Z2t*C1bar';
        [alphahat,betahat,S,rho1]= est_alpha_beta(R0t,R1t,r,Joh_j);
        
        
        % estimate C:
        Chat = -(C0bar(:,1:n) - alphahat*betahat'*C1bar(:,1:n));
        betahat = H(:,1:tH)*betahat; 
       
        % calculate residuals
        res = Z0t  + Z2t*Chat' - Z1t*betahat*alphahat';
         
        switch Joh_j
        %    case 1
        %        res = detrend(res,'linear');
            case 2 % restricted linear trend
                yh=[1:(T-2)]'/(T-2);
        %        res = detrend(res,0);
                res = res - yh*rho1*alphahat';
        %    case 3 % unrestricted constant
        %        res = detrend(res,0);
            case 4 % restricted constant.
                yh=ones(T-2,1);
                res = res - yh*rho1*alphahat';
            otherwise
        end
        
        %LL =  log(det(res'*res/(T-2)));
        
    otherwise % here comes the unrestricted estimation.
        
        % differentiate, if columns are left:
        if r >= tH % no additional columns left
            if (r<tH)
                disp('RH_vecm_H0: Specification is wrong! r<t.');
                return;
            else % no r==t, restriction fully defines cointegrating space -> OLS regression.
                Z1th = Z1t*H(:,1:tH);
                % adjust for deterministics.
                plus1 = 0;
                switch Joh_j
                    case 2
                        yh=[1:(T-2)]'/(T-2);
                        Z1th = [yh,Z1th];
                        plus1 = 1;
                    case 4
                        yh=ones(T-2,1);
                        Z1th = [yh,Z1th];
                        plus1 = 1;
                end;
                
                % calculate regression
                betahat = H(:,1:tH);
                alph = ([Z1th,Z2t]\Z0t);
                alphahat = alph(plus1+[1:tH],:)';
                res = Z0t - [Z1th,Z2t]*alph;
                LL =  (T-2)*log(det(res'*res/(T-2)))+(T-2)*s;
                Chat = -alph(plus1+tH+1:end,:)';
            end
            
            return
        end
        
        Z1th = Z1t*H(:,1:tH);
        
        C0 = (Z2t\Z0t);
        C1 = (Z2t\Z1th);
        
        % generate R_{0t}
        R0t = Z0t - Z2t*C0;
        
        % generate R_{1t}
        R1t = Z1th - Z2t*C1;
        
        [alphahat,betahat,S,rho1]= est_alpha_beta(R0t,R1t,r,Joh_j);
        betahat = H(:,1:tH)*betahat;
        % residuals
        res = Z0t - Z1t*betahat*alphahat';
        Ctilde = Z2t\res;
        Chat = -Ctilde';
        
        % log likelihood
        res = Z0t  + Z2t*Chat' - Z1t*betahat*alphahat';
        
        switch Joh_j
            case 1
                res = detrend(res,'linear');
            case 2 % restricted linear trend
                yh=[1:(T-2)]'/(T-2);
                res = detrend(res,0);
                res = res - yh*rho1*alphahat';
            case 3 % unrestricted constant
                res = detrend(res,0);
            case 4 % restricted constant.
                yh=ones(T-2,1);
                res = res - yh*rho1*alphahat';
            otherwise
        end
        
       % LL =  log(det(res'*res/(T-2)));
end

LL =  (T-2)*log(det(res'*res/(T-2)))+(T-2)*s;



