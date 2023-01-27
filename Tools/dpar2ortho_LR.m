function dC = dpar2ortho_LR(par,n,p,dp);
% delivers a matrix C'C=I using the parameter values in the vector par
%
% SYNTAX: function C = dpar2ortho_LR(par,n,p,dp);
%
% INPUTS:    par ... dx1; parameter vector;
%            n,p ... integers; output matrix C in R^{n x p}.
%            dp  ... integer; used for calculating the derivative with
%            respect to the dp-th entry. 
%
% OUTPUT:    C ... n x p matrix, C = Q_L [I_p;0] Q_R.
%
% REMARK: The matrix is parameterized as Q_L [I;0] Q_R with Q_L and Q_R given
% as products of Givens rotations.
%
% dbauer, 20.11.2019

dC = [eye(p);zeros(n-p,p)];
dp = length(par)-dp+1;
% multiply using the entries of theta_L
for i=1:p
    for j=(n-1):-1:p
        pc = par(end);
        par = par(1:end-1);
        dp = dp-1;
        [Q,in] = dQ(pc,i,j,dp,n);
        %Q= [cos(pc),-sin(pc);sin(pc),cos(pc)];
        dC(in,:)=Q*dC(in,:);
    end;
end;

% multiply from the right using theta_R
for i=(p-1):-1:1
    for j=p:-1:(i+1)
        pc = par(end);
        par = par(1:end-1);
        dp = dp-1;
        [Q,in] = dQ(pc,i,j-1,dp,p);
        %Q= [cos(pc),-sin(pc);sin(pc),cos(pc)]';
        dC(:,in)=dC(:,in)*Q';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% helper function %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q,in] = dQ(pc,i,j,dp,n);

if dp
    in = [i,j+1];
    Q = [cos(pc),-sin(pc);sin(pc),cos(pc)];
else
    in = 1:n;
    Q = zeros(n,n);
    Q([i,j+1],[i,j+1])=[-sin(pc),-cos(pc);cos(pc),-sin(pc)];
end;
