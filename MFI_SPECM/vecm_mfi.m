classdef vecm_mfi
    % MFI model in VAR representation. 
    properties
        urs = zeros(0,3); % unit root structure
        k = 1; % max lag length
        alpha = [];
        beta = [];
        Gamma = [];
        Omega = [];
        s = 1;
    end
    
    methods
        function VAR = conv_vecm_var(vecm_mfi);
            % convert vecm to var matrices.
            s = vecm_mfi.s;
            urs = vecm_mfi.urs;
            M = size(urs,1);
            k = vecm_mfi.k;
            Gammam = vecm_mfi.Gamma;
            for jk=1:k
                Gamma(:,:,jk)=Gammam(:,[(jk-1)*s+1:(jk*s)]);
            end;
            % generate frequencies and p(L)
            p = 1;
            for m=1:M
                pm{m} = 1;
            end;
            
            for m=1:M
                if urs(m,3) % complex root
                    np = [1,-exp(sqrt(-1)*pi*urs(m,1))];
                    p = conv(p,np);
                    p = conv(p,conj(np));
                    for j=1:M
                        if j~=m
                            pm{j} = conv(pm{j},np);
                            pm{j} = conv(pm{j},conj(np));
                        end;
                    end;
                else % real root
                    np = [1,-real(exp(sqrt(-1)*pi*urs(m,1)))];
                    p = conv(p,np);
                    for j=1:M
                        if j~=m
                            pm{j} = conv(pm{j},np);
                        end;
                    end;
                end;
            end;
            K = k+length(p)-1;
            VAR = zeros(s,s,K);
            
            for d=1:length(p)-1 % p(L)X_t
                VAR(:,:,d) = p(d+1)*eye(s); 
            end;
            
            % Gamma_k
            for jk=1:k    
                VAR(:,:,jk) = VAR(:,:,jk) - Gamma(:,:,jk);
                for j=1:length(p)-1
                    VAR(:,:,jk+j) = VAR(:,:,jk+j) - p(j+1)*Gamma(:,:,jk); 
                end;
            end;
            
            % alpha's and beta's.
            for m=1:M
                ah = vecm_mfi.alpha{m};
                bh = vecm_mfi.beta{m};
                if urs(m,3) % complex root
                    fr = (exp(sqrt(-1)*pi*urs(m,1)));
                    pmm = conv(pm{m},[1,-(fr)]);
                    frv = fr.^[1:length(pmm)];
                    cm = urs(m,4)/2;
                    ah_c = (ah(:,1:cm)+sqrt(-1)*ah(:,(cm+1):end))/2;
                    bh_c = (bh(1:s,1:cm)+sqrt(-1)*bh((s+1):end,1:cm));
                    frv = fr.^[1:length(pmm)];
                    
                    AB = ah_c*bh_c'/((pmm*conj(frv)'));
                    for jm=1:length(pmm)
                        VAR(:,:,jm)=VAR(:,:,jm)-2*real(AB*pmm(jm));
                    end;
                    
                else % real root
                    fr = real(exp(sqrt(-1)*pi*urs(m,1)));
                    pmm = pm{m};
                    frv = fr.^[1:length(pmm)];
                    
                    AB = ah*bh'/(pmm*frv');
                    for jm=1:length(pmm)
                        VAR(:,:,jm)=VAR(:,:,jm)-AB*pmm(jm);
                    end;
                end;
            end;
            
            
        end; % function conv_vecm_var.
        
        function IR = cal_IR(vecm_mfi,L);
            % calculate the first L impulse response coefficients.
            
            % tranform into VAR representation
            VAR = conv_vecm_var(vecm_mfi);
            
            % transform to companion state space form
            [s,~,K]=size(VAR);
            A = [zeros(s,s*K);eye(s*(K-1)),zeros(s*(K-1),s)];
            for k=1:K
                A(1:s,(k-1)*s+[1:s])=-VAR(:,:,k);
            end;
            C = A(1:s,:);
            K = [eye(s);zeros(s*(K-1),s)];
                        
            % calculate IR for state space
            IR = zeros(s,s,L);
            for l=1:L 
                IR(:,:,l)= C*A^(l-1)*K;
            end;
            
        end
        
        function plot_IR(vecm_mfi,L);
            % plot the impulse response. 
            
            figure;
            hold on;
            leg = cell(0);
            s = vecm_mfi.s; 
            IR = cal_IR(vecm_mfi,L); 
            
            for a=1:s
                for b=1:s
                    plot(squeeze(IR(a,b,:)));
                    leg{end+1} = sprintf('Row: %d, col: %d',a,b);
                end;
            end;
            legend(leg);
            title('Impulse response function');
        end;
        
        function tf = cal_tf(vecm_mfi,fr);
            % evaluate transfer function at frequencies.
            s = vecm_mfi.s;
            % transform to VAR
            VAR = conv_vecm_var(vecm_mfi);
            [s,~,K]=size(VAR);
            
            tf = zeros(s,s,length(fr));
            for f=1:length(fr)
                om = exp(sqrt(-1)*fr(f)*pi);
                tfh = eye(s);
                for k=1:K
                    tfh = tfh + om^k*squeeze(VAR(:,:,k));
                end;
                tf(:,:,f)=inv(tfh);
            end;
            
        end;
        
        function plot_tf(vecm_mfi,fr);
            % plot the transfer function. 
            
            figure;
            hold on;
            leg = cell(0);
            s = vecm_mfi.s; 
            tf = cal_tf(vecm_mfi,fr);
            
            for a=1:s
                for b=1:s
                    plot(squeeze(abs(tf(a,b,:))));
                    set(gca,'yscale','log');
                    leg{end+1} = sprintf('Row: %d, col: %d',a,b);
                end;
            end;
            legend(leg);
            title('Transfer function');
        end;
    end
end

        