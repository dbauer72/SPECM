function [urs,nums,vals] = test_unit_roots_MFI(Ahat,S,T,s,det_fr);
% test_unit_roots_MFI uses the eigenvalue based test of Bauer and
% Buschmeier to provide estimates of the location of unit roots for the
% estimated matrix Ahat from the CCA procedure.
%
% SYNTAX: [urs,nums,vals] = test_unit_roots_MFI(Ahat,S,T,s,det_fr);
%
% INPUT: Ahat   ... nxn real matrix; describing the dynamics.
%        S      ... integer; number of seasons.
%        T      ... integer; sample size.
%        s      ... integer; dimension of output (acts as an upper bound).
%        det_fr ... Sx1 integer vector; indicators of deterministics at the
%                   S unit root frequencies included in the estimation.
%
% OUTPUT:   urs ... unit root structure; each row corresponds to a unit
% root; [om_j,complex?,c_j]: * om_j is frequency as a fraction of pi.
%                           * complex indicates a complex unit root.
%                           *  c_j denotes the number of common trend.
%           nums ... Sx1 integer vector: c_j numbers for all freqs.
%           vals ... Sx1 cell array containing the absolute distances to
%           the unit root.
%
% REMARK:
% (1) The test compares the sum of the absolute values of the distances
% of the estimated roots to the unit roots with the test statistic as
% described in Bauer and Buschmeier (2021). Real and complex cases are
% separated.
% (2) For each complex unit root, complex conjugated are identical
%     but listed for symmetry reasons.
%
% AUTHOR: dbauer, 7.1.2021.

if length(det_fr)==1
    det_fr = det_fr*ones(S,1);
end;

om = 2*[0:(S-1)]/S; % frequencies.
Sm1 = round(S/2)+1;

nums = zeros(S,1);
vals = cell(S,1);
urs = zeros(0,3);

% calculate eigenvalues
ev = eig(Ahat);

% load critical values
load crit_MFI_ev_test;

for j=1:Sm1
    zj = exp(sqrt(-1)*pi*om(j));
    dist = T*sort(abs(ev-zj),'ascend');
    switch j
        case {1,Sm1} % real unit root
            comp=0;
            if det_fr(j)
                cvals = crit_val_real_det;
            else
                cvals = crit_val_real_wodet;
            end
        otherwise % complex unit root
            comp=1;
            if det_fr(j)
                cvals = crit_val_comp_det;
            else
                cvals = crit_val_comp_wodet;
            end
    end
    cj = seq_trace_test(dist(1:s),cvals(1:s));
    nums(j)=cj;
    if cj>0
        urs(end+1,:)= [om(j),comp,cj];
    end;
    vals{j} = dist(1:s)./([1:s])';
end

for j=(Sm1+1):S % add complex conjugate roots
    nums(j)=nums(S-j+2);
    vals{j}= vals{S-j+2};
end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  seq_trace_test    %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% implements sequential testing: order is increased until the null is not
% rejected any more.
function cj = seq_trace_test(dist,cvals);

i = 1;
dist= cumsum(dist);
s = length(cvals);
while i <= s
    if i == 1
        if dist(s-i+1) <= cvals(s-i+1) % not rejected -> full set of unit roots. 
            cj = s;
            i = s+1;
        end
    end
    if i <= s
        if dist(s-i+1) > cvals(s-i+1)
            i = i + 1;
        else
            cj = s-i+1;
            i = s+1;
        end
        if i == s
            if dist(1) > cvals(1)
                cj = 0;
            else
                cj= 1;
            end
            i = s+1;
        end
    end
end
end
