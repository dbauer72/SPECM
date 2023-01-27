function dZ2t = gen_dZ2t(Z0t,Z2t,Abar,B);
% gen_dZ2t generates the regressors in the C(alpha) procedure for
% estimating alpha. 
%
% SYNTAX: dZ2t = gen_dZ2t(Z0t,Z2t,tA,tB,Chat);
%
% INPUT: Z0t ... Delta y_t concentrated w.r.t. deterministics.
%        Z2t ... tilde x_t concentrated w.r.t. deterministics.
%        Abar  ... Abar 
%        B  ... B
%
% OUTPUT:  dZ2t ... matrix of additional regressors.

[n,s] = size(B);
T = size(Z2t,1); 

Z2tm1 = [ones(1,n);Z2t(1:end-1,:)];
dZ2t = zeros(T,0);

if n<s % few states -> parameters in B and Abar. 
    % nx(s-n) parameters in the last columns of B 
    for a=(s+1):n
        for b=1:n
            dB = zeros(n,s);
            dB(b,a)=1;
            dZ2t = [dZ2t,ltitr(Abar,dB,Z0t)];
        end
    end

    % nxn parameters in Abar
    for a=1:n
        for b=1:n
            dA = zeros(n,n);
            dA(b,a)=1;
            dZ2t = [dZ2t,ltitr(Abar,B,Z2tm1*dA')];
        end
    end
else % more states than dimensions -> parameters only in Abar. 
    for a=1:s
        for b=1:n
            dA = zeros(n,n);
            dA(b,n+1-a)=1;
            dZ2t = [dZ2t,ltitr(Abar,B,Z2tm1*dA')];
        end
    end
end 
