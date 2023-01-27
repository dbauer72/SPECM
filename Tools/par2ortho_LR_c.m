function [C,QL] = par2ortho_LR_c(par,n,p);
% par2ortho_LR_c transforms the parameters in par into the unitary column C.
% Inverse to ortho2par_LR_c. 
%
% SYNTAX: C = par2ortho_LR_c(par,s,c);
%
% INPUT:  par ... dx1 parameter vector.
%         s,c ... dimensions of C. 
%      
% OUTPUT: C ... sxc matrix, C'C=I_c. 
%
% REMARK: calculates C = Q_L [D_c;0]Q_R where Q_L=Q_L(par) and Q_R =
% Q_R(par). 

% AUTHOR: dbauer, 14.10.2020
dd = exp(sqrt(-1)*par(end-p+1:end));
par(end-p+1:end)=[];

C = [diag(dd);zeros(n-p,p)];

% multiply using the entries of theta_L
QL = eye(n);
for i=1:p
    for j=(n-1):-1:p
        Qi = parc2Q(par(end-1:end));
        par(end-1:end)=[];       
        C([i,j+1],:)=Qi'*C([i,j+1],:);
        QL([i,j+1],:)=Qi'*QL([i,j+1],:);        
    end;
end;


% multiply from the right using theta_R
for i=(p-1):-1:1
    for j=p:-1:(i+1)
        Qi = parc2Q(par(end-1:end));
        par(end-1:end)=[];   
        C(:,[i,j])=C(:,[i,j])*conj(Qi);
    end
end

end

function Qi = parc2Q(pa);
% provides the 2x2 Givens complex Givens rotation. 
%
% dbauer, 14.10.2020

i = sqrt(-1);
Qi = [cos(pa(1)),sin(pa(1))*exp(-i*pa(2));-sin(pa(1))*exp(i*pa(2)),cos(pa(1))];
end