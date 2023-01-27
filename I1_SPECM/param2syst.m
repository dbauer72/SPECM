function [A,K,C,D,Omega,th,gr_syst] = param2syst(param,s,n,m,c,restrict,Omega); 
% takes the parameter vector, the integers and the restrictions and
% converts it to system plus innovation variance.
%
% SYNTAX: [A,K,C,D,Omega,th,gr_syst] = param2syst(param,s,n,m,c,restrict); 
%
% INPUT: param ... dx1 parameter vector.
%        s     ... integer; dimension of endogenous vars.
%        n     ... integer; state dimension.
%        m     ... integer; dimension of exogenous vars.
%        c     ... integer; number of common trends.
%        restrict ... structure describing the restrictions. 
%
% OUTPUT: A,K,C,D ... state space system.
%         Omega   ... sxs innovation variance.
%         th      ... thetat structure corresponding to A,K,C,D.
%         gr_syst ... if wanted in the output contains the derivative of
%                       the system matrices with respect to the parameters. 
%                      Array: gr_syst(1).A, .K, .C, .D, .Omega contains the
%                      derivative of A,K,C,D, Omega respectively w.r.t. the
%                      first parameter, gr_syst(2).A, ... w.r.t. to the
%                      second param.
%
% REMARK: uses the parameterization corr. to the param paper. 
%
% AUTHOR: dbauer, 21.11.2019

det_res = restrict.det_res;
if isfield(restrict,'scale')
    scale = diag(restrict.scale); 
else
    scale = eye(length(param));
end;

param = scale*param;

% extract parameters for Omega
if nargin<7 % Omega not contained
    parom = param(1:(s*(s+1)/2));
    nom = length(parom);
    param(1:(s*(s+1)/2))=[];
else
    parom = extr_lowtri(Omega);
end

if nargout>6
    np = length(param)+nom;
    [Omega,dOmega] = fill_lowtri(parom,s);
else
    Omega = fill_lowtri(parom,s);
end

% convert params to matrices
nd = s*m;
if det_res 
    nd = nd- (s-c);
end;

if nargout>6
    nth = length(param);
    [th,param,dth] = param2th(param,s,n,c,restrict);
    nth = nth-length(param);    
else
    [th,param] = param2th(param,s,n,c,restrict);
end

A = th.A;
K = th.K;
C = th.C;
% parameters for deterministics. 
if c>0
    [Q,R]= qr(C(:,1:c));
else
    Q = eye(s);
end;

if (det_res>0)&(c>0)
    D = zeros(s,m);
    pc = param(1:c);
    D(:,1)=C(:,1:c)*pc(:);

    param(1:c)=[];
    if m>1
        D(:,2:end)= reshape(param,s,m-1);
    end
    if nargout>6
        for j=1:c
            dD(:,1,j)=C(:,j);
        end
        if m>1 
            dp = eye(s*(m-1));
            for j=1:(s*(m-1))
                dD(:,:,j) = reshape(dp(:,j),s,m-1);
            end
        end
    end
else
    D = reshape(param,s,m);
    if nargout>6
        dp = eye(s*m);
        for j=1:s*m,
            dD(:,:,j) = reshape(dp(:,j),s,m);
        end
    end;
end

th.A =A;
th.K =K;
th.C =C;
th.D =D;
th.B = zeros(n,m);
th.Omega =Omega;
th.ur = 'I(1)';

if nargout>6 % gradient also wanted    
    for j=1:np
        da = zeros(n,n);
        dk = zeros(n,s);
        dc = zeros(s,n);
        dd = zeros(s,m);
        dom = zeros(s,s);
        
        if j<= nom
            dom = squeeze(dOmega(:,:,j));
        end
        if (j> nom)&(j<=nom+nth) % parameter for (A,K,C)
            da = squeeze(dth.A(:,:,j-nom));
            dc = squeeze(dth.C(:,:,j-nom));
            dk = squeeze(dth.K(:,:,j-nom));
        end
        if j>nom+nth % param for D
            dd = squeeze(dD(:,:,j-nom-nth));
        end        
        
        gr_syst(j).A= da;    
        gr_syst(j).C= dc;    
        gr_syst(j).K= dk;    
        gr_syst(j).D= dd;   
        gr_syst(j).Omega= dom;    
    end
end


