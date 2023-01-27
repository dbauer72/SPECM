function [LL,alphahat,betahat,Chat,V_alpha,V_beta] = RH_specm(y,Abar,B,r,rest,Joh_j,betaquer);
% ribarits_hanzon idea for estimation of RH_SPECM form
%
% SYNTAX: [LL,alphahat,betahat,Chat] = RH_SPECM(y,Abar,B,r,rest,Joh_j);
%
% INPUT: y    ... Txs data matrix,
%        Abar ... nxn matrix A-BC.
%        B    ... nxs matrix
%        r    ... rank of alpha beta' 
%        rest ... 'y'|'n' restricted version.
%        Joh_j ... integer; indicating the version of Johansen's
%        deterministics. 
%        betaquer ... sxr matrix for normalising betahat (without dets).
%
% OUTPUT: LL   ... minimized log likelihood.
%         alphahat, betahat ... sxr matrices.
%         Chat ... sxn optimal C.
%         V_alpha, V_beta ... estimates of variance matrix of alpha and
%         beta respectively. 
%
%
% COMMENT: The restriction C*(I-Abar)^(-1)*B = I+alpha beta' is imposed by
% rest = 'y'.  
%
% AUTHOR: dbauer, 21.11.2017. 

% prepare regressions 
if nargin< 6
    Joh_j = 5; % no deterministic terms. 
end


[T,s] = size(y);
n = size(Abar,1);



if nargin< 7 
    betaquer = zeros(s,r);
    betaquer(1:r,1:r) = eye(r);
end

if size(betaquer,2) ~= r
    betaquer = betaquer(:,1:r);
    disp('Dimension of betaquer needed to be changed!');
end; 

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

% check, if dimensions fit. 
% r >=s-n! Number of states provides a lower bound for the cointegrating rank.  
% find out, if this is possible at all
if (n+r)<s
    disp('r too small, not possible!');
    LL = -Inf;
    alphahat = zeros(s,0);
    betahat = zeros(s+plus1,0);
    Chat = zeros(s,n);
    V_alpha =[];
    V_beta = [];
    return;
end




changedr= 0; 


switch rest
    case 'y'
        % set up restriction
        R = inv(eye(n)-Abar)*B;
        % add restrictions, if not enough
        if n<s % less states than output dimensions
            [LL,alphahat2,betahat2,Chat] = RH_specm(y,Abar,B,r,'n',Joh_j);
            beta2 = betahat2(:,1:(s-n));
            alpha2 = alphahat2(:,1:(s-n));
            R = [R;beta2(1:s,:)'];
            Z2t = [Z2t,Z1t*beta2];
            r = r-s+n;
            disp('Changed r');
            changedr= 1; 
        end;
        if r>0
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
        %        C0bar = C0hat - (C0hat*R+eye(s))*inv(R'*itM*R)*R'*itM;
        %        C1bar = C1hat - (C1hat*R-eye(s))*inv(R'*itM*R)*R'*itM;
        
            % generate R_{0t}
            R0t = Z0t - Z2t*C0bar';
        
            % generate R_{1t}
            R1t = Z1t - Z2t*C1bar';
            [alphahat,betahat,S,rho1]= est_alpha_beta(R0t,R1t,r,Joh_j);
        

            % estimate C:
            if plus1
                Chat = -(C0bar - alphahat*(betahat(plus1+1:end,:)'*C1bar(plus1+1:end,:)+betahat(1,:)'*C1bar(plus1,:)));
            else
                Chat = -(C0bar - alphahat*betahat'*C1bar(plus1+1:end,:));
            end


            %switch Joh_j
            %    case 1
            %        R0t = detrend(R0t,'linear');
            %    case {2,3}
            %        R0t = detrend(R0t,0);
            %end;
        
            % log lokelihood
            %LL = log(det(R0t'*R0t/T))+sum(log(1-(S(1:r).^2)));
            if n<s
                betahat = [betahat,beta2];
                alphahat = [alphahat,alpha2]; 
            end
            % calculate residuals
            res = Z0t  + Z2t*Chat' - Z1t*betahat*alphahat';
            Chat = Chat(:,1:n);
            switch Joh_j
                case 2 % restricted linear trend
        %         res = res - Z1t(:,plus1)*rho1*alphahat';
                
                case 4 % restricted constant.
            %        res = res - Z1t(:,plus1)*rho1*alphahat';
            end                
        else % r=0! 

            if n<s
                alphahat = alpha2;
                betahat = beta2; 
                iR = inv(R);
                Chat = iR(:,1:n);
            elseif n==s
                alphahat = zeros(s,0);
                betahat = zeros(s+plus1,0);
                Chat = inv(R);
            else
                alphahat = zeros(s,0);
                betahat = zeros(s+plus1,0);
                [QR,RR]= qr(R);
                tZ2t = Z2t*QR;
                R0t = Z0t - tZ2t(:,1:s)*inv(RR(1:s,:))';
                R2t = tZ2t(:,(s+1):end);
                Chat1 = (R2t\R0t)'; 

                Chat = inv(RR(1:s,1:s))*QR(1:n,1:s)'+Chat1*QR(:,(s+1):end)';
            end
            res = Z0t - Z2t(:,1:n)*Chat'; 
        end
    otherwise % here comes the unrestricted estimation.
        C0 = (Z2t\Z0t);
        C1 = (Z2t\Z1t);
        
        % generate R_{0t}
        R0t = Z0t - Z2t*C0;
        
        % generate R_{1t}
        R1t = Z1t - Z2t*C1;
        
        [alphahat,betahat,S,rho1]= est_alpha_beta(R0t,R1t,r,Joh_j);
        
        %switch Joh_j
        %    case 1
        %        R0t = detrend(R0t,'linear');
        %    case {2,3}
        %        R0t = detrend(R0t,0);
        %end;
        % residuals
        
        
        % log likelihood
        if r>0
            res = Z0t  - Z1t*betahat*alphahat';
        else
            res = Z0t;
        end; 

        switch Joh_j
            case 2 % restricted linear trend
               % res = res - Z1t(:,plus1)*rho1*alphahat';                
            case 4 % restricted constant.
                %res = res - Z1t(:,plus1)*rho1*alphahat';
        end
        
        Ctilde = Z2t\res;
        Chat = -Ctilde';
        
        res = res - Z2t*Ctilde;
end

%% 
[mm,nn]=size(betahat);
try 
    bb = betaquer.*betahat;
catch
    betaquer = [eye(nn);zeros(mm-nn,nn)];
end

Trafo = betaquer'*betahat; 
betahat = betahat*inv(Trafo);
alphahat = alphahat * Trafo';

Omegahat = res'*res/(T-2);
if r>0 && (changedr ==0)
    Sig_bb = betahat'*R1t'*R1t*betahat; 

    % Johansen uses Kronecker products for rowwise vectorization!
    % V_alpha = kron (Sigma_bb^(-1), Omega), see p. Johansen p. 182. 
    V_alpha = kron(inv(Sig_bb),Omegahat);

% V_beta = kron( inv(int GG')', (alpha' Omega^(-1) alpha)^{-1}), see
% Johansen p. 181.


    S11 = R1t'*R1t; 

    Proj = eye(mm) - betahat*betaquer'; 
    V_beta = kron(inv(alphahat'*inv(Omegahat)*alphahat),Proj*inv(S11)*Proj'); 
else
    V_beta = zeros(0,0);
    V_alpha = zeros(0,0);
end;
LL =  (T-2)*log(det(res'*res/(T-2)))+(T-2)*s;




