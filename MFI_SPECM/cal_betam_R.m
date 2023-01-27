function [beta,LL] = cal_betam_R(S00,S01,S11,r);
% Estimates the parameters in the partial loglikelihood for one root only. 
%
% see Johansen, Schaumburg, 1999, p. 319, (16).
% REMARK: implements the restriction of real valued beta. 
% AUTHOR: Martin Wagner.

s = size(S00,1);

S11d0 = S11 - (S01'/S00)*S01;

% rotate such that heading beta with I is not restrictive. 
sn = size(S11,1)/2;
[v,d] = seig(S11(1:sn,1:sn));
TT = [v(1:sn,sn:-1:1),zeros(sn,sn);zeros(sn,sn),v(1:sn,sn:-1:1)];

tS11 = TT'*S11*TT;
tS01 = S01*TT;
tS11d0 = TT'*S11d0*TT;

ab0 = tS01*inv(tS11); % unrestricted estimate
[U,S,V]= svd(ab0);

if (r<s) &&(s>1) && (r>0)
    
    % convert to complex for starting value 
    TT2 = [eye(sn),sqrt(-1)*eye(sn)];
    iS11 = TT2*tS11*TT2';
    iS01 = tS01*TT2';
    icc = iS01*inv(iS11);
    % svd for rrr in complex terms
    [U,S,V]= svd(icc);
    ibb = V(:,1:r) * inv(V(1:r,1:r));
    % convert back to real terms
    bb1n = [real(ibb(r+1:end,:))];

    % numerical optimization
    [betamn,LL] = fminunc(@(x) ll_det_R(x,tS11,tS11d0,r),bb1n(:));
    % convert back. 
    beta = TT*cal_beta_par_R(sn,r,betamn);
    LL = (LL + log(det(S00)));
else
    if (r==s)
        beta = eye(2*r);
        LL = (log(det(S00))+log(det(S11d0))-log(det(S11)));
    end;
end;

%%
function dd = ll_det_R(x,S11,S11d0,r);
% logdet det(b'S11d0 b)/det(b'S11 b).

s = size(S11,1)/2;
beta = cal_beta_par_R(s,r,x);
dd = log(det(beta'*S11d0*beta)/det(beta'*S11*beta));

%%
function Betam = cal_beta_par_R(s,r,x);
%
bb1 = reshape(x,s-r,r);
bb0 = [eye(r);bb1(:,1:r)];
Betam = [real(bb0),-imag(bb0);imag(bb0),real(bb0)];


  