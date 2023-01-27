function th = VAR(y,u,p,which);
%
% 	estimate a VAR(X) in levels with lag length p.
%
% SYNTAX:  th = VAR(y,u,p);
%
% INPUTS:   y ... Tx s; observations of dependent var (y) 
%           u ... T x m; observations of exogenous var (u)
%		    p ... lag length
%           which ... 'LS'|'YW' 
%
% OUTPUTS:  
%          th ... theta structure corresponding to VARX(kp).
%		  
% REMARK: 'LS' uses least squares estimation, 'YW' is based on the Yule
% Walker equations. 
%
%  AUTHOR: adaptations by dbauer, 18.8.2020.

[T,s] = size(y);
nu = size(u,2);

if nargin<4
    which = 'LS';
end



if which == 'LS' % LS estimation
    Z = u(p+1:end,:);
    for j=1:p
        Z = [Z,y(p-j+[1:(T-p)],:)];
    end;
    phi = (Z\y(p+1:end,:))';
    res = y(p+1:end,:) - Z*phi';
    Omega = res'*res/(T-p);
else % Yule Walker equations
    z = [y,u];
    nz = s+nu;
    R = mcovf(z,p+1); % mcovf... estimates covariance sequence.
    R=reshape(R,nz,(p+1)*nz);
    
    Gm=zeros(p*s+nu,p*s+nu);
    H=R(1:s,s+[1:nu]);
    Gm(1:nu,1:nu) = R(s+[1:nu],s+[1:nu]);
    
    for i=1:p
        H = [H,R(1:s,i*nz+[1:s])];
        Gm(1:nu,nu+(i-1)*s+[1:s]) = R(s+[1:nu],i*nz+[1:s]);
        Gm(nu+(i-1)*s+[1:s],1:nu) = R(s+[1:nu],i*nz+[1:s])';
        for j=i:p
            Gm(nu+(i-1)*s+[1:s],nu+(j-1)*s+[1:s])= R(1:s,(j-i)*nz+[1:s]);
            Gm(nu+(j-1)*s+[1:s],nu+(i-1)*s+[1:s])= R(1:s,(j-i)*nz+[1:s])';
        end;
    end;
    phi = H*inv(Gm);
    Omega = R(1:s,1:s)-phi*H';
end;

% fill in matrices into theta format. 
th = theta_urs();
th.which = 'poly'; %use polynomial form.
th.a = [eye(s),-phi(:,nu+1:end)];
th.d = phi(:,1:nu);
th.b = eye(s);
th.m = nu;

th.Omega = Omega;



    
