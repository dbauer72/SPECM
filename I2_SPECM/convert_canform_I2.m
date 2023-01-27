function [A,K,C] = convert_canform_I2(ae,ke,ce,ds); 
% converts a general state space system into the Bauer Wagner canonical
% form in the I(2) case. 
%
% SYNTAX: [Ae,Ke,Ce] = convert_canform_I2(ae,ke,ce,ds); 
%
% INPUT: ae,ke,ce ... state space system
%        ds        ... 2x1 pair of integers, ds(1): number of I(2) common
%                      trends, ds(2) ... number of additional I(1) common trends. 
%
% OUTPUT: Ae,Ke,Ce ... state space system
%
% AUTHOR: dbauer, 13.2.2020.

% integers 
[n,s] = size(ke);
c = ds(1)*2+ds(2);

nbull = n-c; 

% start with Jordan normal form
[Tv,A] = Jordan_Au(ae,ds);

%A = Tv*ae*inv(Tv);
K = Tv*ke;
C = ce*inv(Tv);

% impose other restrictions in canonical form
% start with C1.
C1 = C(:,1:ds(1));
[q,r]=qr(C1);
C(:,1:ds(1)) = q(:,1:ds(1));
K(1:ds(1),:) = r(1:ds(1),1:ds(1))*K(1:ds(1),:);
K(ds(1)+[1:ds(1)],:) = r(1:ds(1),1:ds(1))*K(ds(1)+[1:ds(1)],:);
C(:,ds(1)+[1:ds(1)]) = C(:,ds(1)+[1:ds(1)])*inv(r(1:ds(1),1:ds(1)));
% next is C3.
Tv = eye(n);
Tv(1:ds(1),2*ds(1)+[1:ds(2)])=-C(:,1:ds(1))'*C(:,2*ds(1)+[1:ds(2)]);
iTv = inv(Tv);
A = iTv*A*Tv;
K = iTv*K;
C = C*Tv;
C3 = C(:,2*ds(1)+[1:ds(2)]);
[q,r]=qr(C3);
C(:,2*ds(1)+[1:ds(2)]) = q(:,1:ds(2));
K(2*ds(1)+[1:ds(2)],:) = r(1:ds(2),1:ds(2))*K(2*ds(1)+[1:ds(2)],:);

% finally C2.
Tv = eye(n);
Tv(2*ds(1)+[1:ds(2)],ds(1)+[1:ds(1)])=-C(:,2*ds(1)+[1:ds(2)])'*C(:,ds(1)+[1:ds(1)]);
Tv([1:ds(1)],ds(1)+[1:ds(1)])=-C(:,[1:ds(1)])'*C(:,ds(1)+[1:ds(1)]);
iTv = inv(Tv);
A = iTv*A*Tv;
K = iTv*K;
C = C*Tv;

% now standardize integrated part
K2 = K(ds(1)+[1:ds(1)],:);
[Q2,R] = qr(K2);

K3 = K(2*ds(1)+[1:ds(2)],:);
[Q3,R] = qr(K3);

Tv =eye(n);
Tv([1:ds(1)],[1:ds(1)])=Q2;
Tv(ds(1)+[1:ds(1)],ds(1)+[1:ds(1)])=Q2;
Tv(2*ds(1)+[1:ds(2)],2*ds(1)+[1:ds(2)])=Q3;
iTv = inv(Tv);
A = iTv*A*Tv;
K = iTv*K;
C = C*Tv;

% standardize stable part
Abull = A(c+1:end,c+1:end);
Kbull = K(c+1:end,:);
Cbull = C(:,c+1:end);

[~,Ab,Cb,Kb] = ss2ech_n(Abull',Cbull',Kbull');

% fill in

A(c+1:end,c+1:end)=Ab';
K(c+1:end,:) = Kb';
C(:,c+1:end) = Cb';
