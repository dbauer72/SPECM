function [tildeP,p,z] = test_season_port(y,Sh,M)
% test_season_port implements the Portmanteau test for the existence of
% seasonal unit roots of Buschmeier (2022) extending the portmanteau test
% of Zhang and Chan. 
%
% SYNTAX: [tildeP,p] = test_season_port(y,S)
%
% INPUT:
%         y ... T x p real matrix of time series
%         Sh ... integer; number of distinct unit roots. Default: Sh=3. 
%         M ... integer; number of covariances to use. Default: estimated
%                using AIC. 
%
% OUTPUT: 
%         tildeP ... Shxp ; real matrix; test results for each unit root (S) and each time
%                    series (p); real matrix; p-values for test statistics.
% 
%         p      ... Sh x p; real matrix; p-values for tests. 
%         z      ... Sh x 2; matrix of unit root frequencies. 
%
% AUTHOR: dbauer, 17.1.2023.

load test_port_cdf 

if nargin<3
    M =100;
end
if nargin<2
    Sh = 3;
end

[T,s]= size(y);

% number of seasons 
S = 2*(Sh-1);
Ts = floor(T/S);
% calculate roots  
z = zeros(Sh,1);
for j=1:Sh
    z(j,1)= 2*pi*(j-1)/S; % frequency
    if (j>1)&&(j<Sh) 
        z(j,2) = 1; % complex unit root
    end
end

tildeP = zeros(Sh,s);
p=tildeP; 

for j=1:Sh  %cycle over seasonal unit roots 
    % filter out all other roots 
    y_f = y; 
    for jf = 1:Sh
        if (jf ~= j) 
            if z(jf,2)==1
                y_f = filter([1,-2*cos(z(jf,1)),1],1,y_f,[],1);
            else
                if (jf==1)
                    y_f = filter([1,-1],1,y_f,[],1);
                else
                    y_f = filter([1,1],1,y_f,[],1);
                end
            end
        end
    end

    for a = 1:s % cycle over coordinates
        y_cur = y_f(:,a);

        if z(j,2)==1 % for complex roots -> normalize. 
            % calculate the vector of seasons repr to estimate long run variance
            Y_tau = reshape(y_cur(1:Ts*S),S,Ts)';

            % estimate long run variance matrix. 
            dY_tau = detrend(Y_tau(2:end,:)-Y_tau(1:end-1,:),1);
            [k,~,~,~,~,~,thar] = aicest(dY_tau,S,-M);
            A = thar.a;
            Omega = thar.Omega; 
            Pi = eye(S);
            for jk=1:k
                Pi = Pi + A(:,jk*S+(1:S));
            end
            Sigmainf = inv(Pi)*Omega*inv(Pi)'; 
            
            % SVD to obtain first two eigenvalues. 
            [U,Lam] = svd(Sigmainf);
            Lamdagger = eye(S);
            Lamdagger(1,1)=1/sqrt(Lam(1,1));
            Lamdagger(2,2)=1/sqrt(Lam(2,2));
            
            % normalize by Lambda_c
            Trans = U* Lamdagger * U';
            Y_tau = Y_tau * Trans; 
            y_cur = reshape(Y_tau',Ts*S,1); 
        end

        % calculate test statistic 
        R=mcovf(y_cur,M); % covariance sequence
        CO = R(1,1:S:end)/R(1,1); % correlation sequence
        MS= length(CO)-1;
        tildeP(j,a) = Ts*sum(1-CO)/(MS*(MS+1)/2);

        % find p-values 
        pp = sum(tildeP(j,a)<test_cdf)/length(test_cdf);
        p(j,a) = pp;
    end
end



end % function 