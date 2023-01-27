function [Ae,Ke,Ce] = convert_canform_I1(ae,ke,ce,c); 
% converts a general state space system into the Bauer Wagner canonical
% form in the I(1) case. 
%
% SYNTAX: [Ae,Ke,Ce] = convert_canform_I1(ae,ke,ce,c); 
%
% INPUT: ae,ke,ce ... state space system
%        c        ... integer, number of common trends.
%
% OUTPUT: Ae,Ke,Ce ... state space system
%
% AUTHOR: dbauer, 28.11.2017.

[d,t]=seig(ae);

% transform to diagonal A
Ae = inv(d)*ae*d;
Ke = inv(d)*ke;
Ce = ce*d;

% parts corresponding to unit roots 
Ae(1:c,1:c)=eye(c);
Ae(1:c,c+1:end)=0;
Ae(c+1:end,1:c)=0;

% convert Ke and Ce. 
CK = real(Ce(:,1:c)*Ke(1:c,:));
[u,s,v]=svd(CK);

Ce(:,1:c) = u(:,1:c);
Ke(1:c,:) = s(1:c,1:c)*v(:,1:c)';
[q,r]=qr(Ke(1:c,1:c));
Ce(:,1:c)=Ce(:,1:c)*q; 
Ke(1:c,:) = q'*Ke(1:c,:);

for j=1:c
    sk = sign(Ke(j,j));
    Ke(j,:)=sk*Ke(j,:);
    Ce(:,j)=Ce(:,j)*sk;
end;

% stable part of transfer function. 
% build impulse response sequence. 
[s,n]=size(Ce(:,c+1:end));

Abullet = Ae(c+1:end,c+1:end);
Cbullet = Ce(:,c+1:end);
Kbullet = Ke(c+1:end,:);

[Ab,Cb,Kb]= ss2ech_n(Abullet,Cbullet,Kbullet);

Ae(c+1:end,c+1:end) = Ab;
Ke(c+1:end,:)= Kb;
Ce(:,c+1:end) = Cb; 

