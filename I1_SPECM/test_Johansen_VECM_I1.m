function [LL,b_eta,alphahat,betahat,th] = test_Johansen_VECM_I1(tily,k,alpha,Joh_j,plots);
% Calculates the Johansen test using the VECM representation 
%
% SYNTAX: [LL,b_eta,alphahat,betahat,Phi] = test_Johansen_VECM(y,k,alpha,Joh_j,plots);
%
% INPUT: y  ... Txs observations
%        k  ... integer; lag length (if negative estimated).
%        alpha ... real; significance level.
%        Joh_j ... integer; 1-5 deterministics according to Johansen.
%        plots ... indicator, if plots are requested. 
%
% OUTPUT: LL      ... (s+1)x1 vector of test results H(j) vs. H(s).
%         b_eta   ... integer; estimated cointegrating rank.
%         alphahat ... sx b_eta matrix of loadings.
%         betahat  ... sx b_eta cointegrating vectors.
%         th       ... theta_structure for estimate using b_eta. 
%
% AUTHOR: dbauer, 19.12.2019

s = size(tily,2); 

for r=0:s,
    [LL(r+1),alphahat,betahat,res,Phi] = VECM_I1(tily,k,r,Joh_j);
end;

% Table really corresponding to model without constant
load eta_critical % MAT file contains simulated critical values. 
etamat = eta_est{Joh_j}; 
% [248.771 256.231 262.694 270.466 
%          207.212 214.069 220.000 226.945
%		  169.302 175.472 181.199 187.821
%          135.164 140.736 145.797 151.696 
%          104.933 109.932 114.241 119.578 			   
%		  78.303 82.614 86.363 91.115 
%          55.54 59.235 62.68 66.705
%		  36.579 39.71 42.588 45.995 
%          21.581 24.076 26.422 29.194
%		  10.347 12.212 13.939 16.158
%          2.98 4.141 5.302 7.015];
      
% Selection of column corresponding to ALPHA
if alpha == 0.10
    col = 1;
elseif alpha == 0.05
    col = 2;
elseif alpha == 0.025
    col = 3;
elseif alpha == 0.01
    col = 4;
end;

% Testing sequence
% TRACE TEST
T = size(tily,1);
tracevec = (LL-LL(end));

i = 1;
while i <= s,
    if i == 1,
        if tracevec(i) <= etamat(11-s+i,col),
            b_eta = 0;
            i = s+1;
        end;
    end;
    if i <= s,
    if tracevec(i) > etamat(11-s+i,col),
        i = i + 1;
    else 
        b_eta = i-1;
        i = s+1;
    end;    
    if i == s,
        if tracevec(i) > etamat(11-s+i,col),
            b_eta = s;
            i = s+1;
        end;            
    end;    
    end;
end;



if plots
    figure;
    hold on;
    set(gca,'XTick',0:s-1);
    bar(0:s-1,etamat(11-s+1:11,col),'b');
    plot(0:s-1,tracevec(1:end-1),'r*');
    title(sprintf('Estimated coint rank: %d',b_eta));
end;

% perform estimation again with c= s-b_eta. 
[~,alphahat,betahat,~,~,th] = VECM_I1(tily,k,b_eta,Joh_j);


