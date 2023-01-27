function [Y,A,B,C] = simu_seas_ur(gamma,rho,sig,T);
%
% model of Cubadda and Omzigt.
%

Y=zeros(2,T+100);
E = randn(T+100,2)';
Sigma = [1,sig*rho;sig*rho,sig^2 ];
C = chol(Sigma)';
E = C*E;

A1 = [-.2;0]*[1,-1];
A2 = [.2;0]*[1,-1];
A3 = [gamma;0]*[1,0];
A4 = [gamma;0]*[0,-1];

% for t=5:(T+100)
%     Y(:,t)=E(:,t)+Y(:,t-4)+A1*(Y(:,t-1)+Y(:,t-2)+Y(:,t-3)+Y(:,t-4))+A2*(Y(:,t-1)-Y(:,t-2)+Y(:,t-3)-Y(:,t-4))+A3*(Y(:,t-1)-Y(:,t-3))+A4*(Y(:,t-2)-Y(:,t-4));    
% end;


A1 = [gamma,0;0,0];
A2 = [-0.4,-gamma+.4;0,0];
A3 = [-gamma,0;0,0];
A4 = [0.6-gamma/10,.4+gamma;0,1];

Y=zeros(2,T+100);

for t=5:(T+100)
    Y(:,t) = A1*Y(:,t-1)+A2*Y(:,t-2)+A3*Y(:,t-3)+A4*Y(:,t-4)+E(:,t);
end;

Y = Y(:,101:end);

A = [A1,A2,A3,A4;eye(6),zeros(6,2)];
C = A(1:2,:);
B = [eye(2);zeros(6,2)];
eig(A)