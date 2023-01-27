function [Tv,Av] = Jordan_Au(A,ds);
% calculates the rotation to transform A into reordered Jordan normal form
% according to I(2) canform.
%
% SYNTAX: [Tv,Av] = Jordan_Au(A,ds);
%
% INPUT: A ... nxn matrix 
%        ds ... 2x1 integer indicating sizes of blocks.
%      
% OUTPUT: Tv ... nxn real matrix of transformation.
%         Av ... nxn real matrix; Av = Tv*A*inv(Tv)
%
% AUTHOR: dbauer, 7.2.2020

n = size(A,1); 
c = 2*ds(1)+ds(2);
nbull = n-c;

% 1) separate unit roots from stationary components
AmI = A - eye(n);
AmI_sq = AmI^2;
[V,D]= eig(AmI_sq);
d = diag(D);
[~,I]=sort(abs(d));



V = V(:,I);
% --- find and eliminate complex values ---
abs_im = sum(abs(imag(V)));
for j=1:n
    if abs_im(j)>10^(-10) % complex valued column found 
        % find second column of same 
        [~,in] = min(sum(abs(V-conj(V(:,j)))));
        V(:,[j,in])= V(:,[j,in])*[1,1;-i,i]';
        abs_im = sum(abs(imag(V)));
    end
end

Tv = inv(V);
Av = Tv*A*inv(Tv);

% 2) deal with unit root part.
Au = Av(1:c,1:c);
AumI = Au - eye(c); % ker(Au -I)= S1 union S3, its image S2
[Uu,Su,Vu]=svd(AumI);  

S2 = Vu(:,1:ds(1));
S1 = AumI*S2; 
S12 = [S1,S2];

% find out, if these two subspaces are really similar. 
si = norm(S2 - S1*(S1\S2));
if si<0.0001 % both are very similar. 
    rr = eye(size(S1,1)) - S1*inv(S1'*S1)*S1';
    S2 = rr(:,1:ds(1));
end
S3t = eye(c) - S12*inv(S12'*S12)*S12';
[q,~,~]=svd(S3t);
S3 = q(:,1:ds(2));

Tvu = inv([S1,S2,S3]);

Tv = [Tvu,zeros(c,nbull);zeros(nbull,c),eye(nbull)]*Tv;
Av = Tv*A*inv(Tv);

% impose structure corresponding to I(2) processes.
Av(1:ds(1),1:ds(1)) = eye(ds(1));
Av(ds(1)+[1:sum(ds)],ds(1)+[1:sum(ds)]) = eye(sum(ds));
Av(ds(1)+[1:sum(ds)],2*ds(1)+ds(2)+1:end) =0;
Av(1:ds(1),ds(1)+[1:ds(1)])=eye(ds(1));
Av(1:ds(1),(2*ds(1)+1):end)=0;
Av(ds(1)+1:end,1:ds(1))=0;
Av(2*ds(1)+1:end,ds(1)+[1:ds(1)])=0;
Av(2*ds(1)+ds(2)+1:end,2*ds(1)+[1:ds(2)])=0;



