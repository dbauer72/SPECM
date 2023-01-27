function [pXk,Xm2] = cal_reg_det(S,pXk,Xm2,urs_full,Joh_j);
% cal_reg_det provides the enhancement due to deterministics for MFI(1)
% Johansen Schaumburg procs. 
%
% SYNTAX: [pXk,Xm2,kwd] = cal_reg_det(S,pXk,Xm2,urs_full,Joh_j);
% 
% INPUT: S          ... integer; number of observations per year
%        pXk        ...   Tx m real; stationary regressors.
%        Xm2        ... T x m x k real; non-stationary filtered terms X_t^{(m)}
%        urs_full   ... unit root structure enanced with all frequencies.
%        Joh_j      ... Johansen category for deterministics. 
% 
% OUTPUT:  pXk ... adjusted matrix
%          Xm2 ... adjusted regressors. 
%          kwd ... number of regressors not integrated.
%
% AUTHOR: dbauer, 20.8.2020.

om = [0:(S/2)]'/S*2; % frequencies between [0,1]*pi.
fr = exp(om*pi*sqrt(-1)); % corr. unit roots. 

s = round(size(Xm2,2)/2);

T = size(pXk,1);

M = size(urs_full,1);

% for m=1 do it differently 
switch Joh_j
    case {1,2,3} % unrestricted constant included 
        pXk = [pXk,ones(T,1)];
    case 4 % restricted constant added to non-stationary regressors.
        Xm2(:,2*s+1,1)=1;    
end

for m=2:M-1 % cycle over unit roots 
    if Joh_j<5 % only for 5 not common cycles
        zm = fr(m).^([1:T]');
        xx = squeeze(Xm2(:,:,m));
        Xm2(:,2*s+2,m)=0;
        Xm2(:,:,m)=[xx(:,1:s),real(zm),xx(:,s+[1:s]),imag(zm)];
%        Xm2(:,2*s+2,m)=imag(zm);        
    end
end

% for z=-1
if Joh_j < 5
    Xm2(:,2*s+1,M)= (-1).^([1:T]');
end
