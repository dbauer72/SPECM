function [LR,b_eta] = test_Johansen_SPECM_MFI1(tily,Abar,K,alpha,S,Joh_j,plots);
% Calculates the Johansen test using the SPECM representation 
%
% SYNTAX: [LL,b_eta,alphahat,betahat,Phi] = test_Johansen_VECM(y,Abar,K,alpha,freq,S,Joh_j,plots);
%
% INPUT: y  ... Txs observations
%        k  ... integer; lag length (if negative estimated).
%        alpha ... real; significance level.
%        S     ... integer; number of observations per year.
%        Joh_j ... integer; 1-5 deterministics according to Johansen (see
%                  VECM_MFI1 for details)
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

% provide unit roots: 
urs(:,1)=[0:(2/S):1]; 
urs(:,2)=1;

fr = urs(:,1);
for j=1:length(fr) 
    if (fr(j)<0.00001)||(fr(j)>0.99999) % frequency very close to 0 or 1 -> real unit roots.  
        urs(j,2)=0;
    end;    
end
urs(:,3)=1; % first specification for estimation: one common trend everywhere. 

% start the estimation. This provides likelihood ratios for all unit roots
% in the cell array LR. 
[LR] = RH_SPECM_MFI1(tily,S,Abar,K,Joh_j);

% Table really corresponding to model without constant
%
load eta_critical_MFI % MAT file contains simulated critical values. 
etamat_real = eta_est{Joh_j};
etamat_complex = eta_est_com{Joh_j};

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% just a replacement, needs to be updated with real values later on!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

if plots 
    figure;
end;

% now iterate over frequencies.
M = length(fr);

for m = 1:M
    if (m==1)||(m==M)
        etamat = etamat_real;
    else
        etamat = etamat_complex;
    end;
    LL = LR{m};
    tracevec = (LL(end)-LL);
    
    i = 1;
    while i <= s,
        if i == 1,
            if tracevec(i) <= etamat(11-s+i,col),
                b_eta(m) = 0;
                i = s+1;
            end;
        end;
        if i <= s,
            if tracevec(i) > etamat(11-s+i,col),
                i = i + 1;
            else
                b_eta(m) = i-1;
                i = s+1;
            end;
            if i == s,
                if tracevec(i) > etamat(11-s+i,col),
                    b_eta(m) = s;
                    i = s+1;
                end;
            end;
        end;
    end;
    if plots
        subplot(M,1,m);
        hold on;
        set(gca,'XTick',0:s-1);
        bar(0:s-1,etamat(11-s+1:11,col),'b');
        plot(0:s-1,tracevec(1:end-1),'r*');
        title(sprintf('Estimated coint rank at frequency %1.2f pi: %d',fr(m),b_eta(m)));
    end;    
end



