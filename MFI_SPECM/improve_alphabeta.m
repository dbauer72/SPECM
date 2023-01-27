function [estimate,LR] = improve_alphabeta(R0t,R1t,estimate,m,which,Joh_j);
% seperates out all but one unit root in the VECM representation of
% Johansen and Schaumburg based on the estimate. Different ideas are
% distinguished according to which.
%
% SYNTAX: [U0t,U1t] = sep_all_but_one(R0t,R1t,estimate,m,which);
%
% INPUT: R0t ... T x s     ... real observations p(l)X_t
%        R1t ... T x s x M ... regressors X_t^(m)
%        estimate ... current estimate.
%        m        ... unit root to focus on.
%        which    ... 'project'|'subtract'
%
% OUTPUT: estimate ... updated estimate.
%         LR       ... vector of JS trace tests. 
%
% AUTHOR: dbauer, 26.11.2015.

[T,s] = size(R0t);


urs = estimate.urs; % unit root structure
cm = s-urs(m,3); % 3rd column: number of common cycles.
tm = urs(m,4);
in = urs(:,4);

M = size(urs,1);



switch which
    case 'project' % project out the other regressors.
        X = zeros(T,0);
        
        for j=1:M
            if j~= m
                if urs(j,2) % complex unit root
                    X = [X,squeeze(R1t(:,:,j))];
                else
                    X = [X,squeeze(R1t(:,1:s,j))];
                    if urs(j,1)<0.00001
                        if Joh_j == 4
                            X(:,end+1)=1;
                        end
                    else
                        if Joh_j<5 
                            X = [X,squeeze(R1t(:,2*s+1,j))];
                        end;
                    end                        
                end;
            end;
        end;
        
        U0t = R0t - X * (X\R0t);
        U1t = squeeze(R1t(:,:,m)) - X * (X\squeeze(R1t(:,:,m)));
        
        
    otherwise
        U1t = squeeze(R1t(:,:,m));
        U0t = R0t;
        
        for j=1:M
            if j~=m
                if urs(j,2) % complex unit root
                    X = squeeze(R1t(:,:,j));
                    Pi = estimate.alpha{urs(j,4)}*estimate.beta{urs(j,4)}';
                    U0t = U0t - X*Pi';
                else
                    X = squeeze(R1t(:,1:s,j));
                    if urs(j,1)<0.00001
                        if Joh_j == 4
                            X(:,end+1)=1;
                        end
                    else
                        if Joh_j<5 
                            X = [X,squeeze(R1t(:,2*s+1,j))];
                        end;
                    end                        

                    Pi = estimate.alpha{urs(j,4)}*estimate.beta{urs(j,4)}';
                    U0t = U0t - X*Pi';
                end;
            end;
        end;
end

S00 = U0t'*U0t/T;
S01 = U0t'*U1t/T;
S11 = U1t'*U1t/T;

switch urs(m,2)
    case 1 % complex root
        me = min(eig(S11));
        if me<10^(-4)
            S11 = S11+eye(size(S11,1))*(10^(-5)-me);
        end;
        me = min(eig(S00));
        if me<10^(-4)
            S00 = S00+eye(s)*(10^(-5)-me);
        end;
        S11d0 = S11 - S01'*inv(S00)*S01;
        me = min(eig(S11d0));
        if me<10^(-6)
            S11d0 = S11d0+eye(size(S11,1))*(10^(-5)-me);
        end;
        if cm>0
            beta = cal_betam(S00,S01,S11,cm);
            U1tb = U1t*beta;
            alpha = (U1tb\U0t)';
            estimate.alpha{tm} = alpha;
            estimate.beta{tm} = beta;
        else
            estimate.alpha{tm} = estimate.alpha{tm}*0;
        end;
        
        % calc LR for test
        LR(1) = real(-T*(log(det(S00))));

        for j=1:(s-1)
            Betam = cal_betam(S00,S01,S11,j);
            LR(j+1) =  real(-T*(log(det(S00))+log(det(Betam'*S11d0*Betam)/det(Betam'*S11*Betam))));
        end;
        LR(s+1) = real(-T*(log(det(S00))+log(det(S11d0)/det(S11))));
        
    otherwise % real root
        %reduce dim
        ind = 1:s; 
        if (urs(m,1)<0.0001)
            if (Joh_j == 4)
                ind = [ind,2*s+1];
            end
        else
            if Joh_j<5
                ind = [ind,2*s+1];
            end
        end
        
        S01 = S01(:,ind);
        S11 = S11(ind,ind);
        S11d0 = S11 - S01'*(S00(1:s,1:s)\S01);
        S11d0(isnan(S11d0))=0;
        S11d0(isinf(S11d0))=0;
        me = min(eig(S11));
        if me<10^(-4)
            S11 = S11+eye(size(S11,1))*(10^(-4)-me);
            S11d0 = S11d0+eye(size(S11,1))*(10^(-4)-me);
        end;
        me = min(eig(S00));
        if me<10^(-4)
            S00 = S00+eye(size(S00,1))*(10^(-4)-me);
        end;
        me = min(eig(S11d0));
        if me<10^(-4)
            S11d0 = S11d0+eye(size(S11,1))*(10^(-4)-me);
        end;
        
        
        c11 = chol(S11);
        c00= chol(S00);
        B01 = inv(c00)'*S01*inv(c11);
        [U,S,V] = svd(B01);
        alphah = c00'*U(:,1:cm)*S(1:cm,1:cm);
        betah = inv(c11)*V(:,1:cm);
        % renormalize
        BB = betah(1:cm,:);
        betah = betah*inv(BB);
        alphah = alphah*BB';
        estimate.alpha{tm} = alphah;
        estimate.beta{tm} = betah;
        
        LR(1) = real(-T*(log(det(S00))));
        for j=1:s
            Betam = inv(c11)*V(:,1:j);
            LR(j+1) = real(-T*(log(det(S00))+log(det(Betam'*S11d0*Betam)/det(Betam'*S11*Betam))));
        end;

end
