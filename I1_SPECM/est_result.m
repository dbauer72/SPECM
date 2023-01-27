classdef est_result
    % results from estimating coint- state space systems.
    %
    % Fields:
    %             theta:        theta structure containing the estimated system
    %             deviance:     -2/T log Gaussian likelihood.
    %             res:          Txs matrix of residuals.
    %             param:        dx1 parameter vector
    %             call:         string; info on call for estimation.
    %             restrict:     structure; info on restrictions.
    %             Pbull:        indicator; if Pbull>0: likelihood is
    %                           calculated for stationary distribution of initial state.
    %             urs:          unit root structure. Implementation depends
    %                           on case.
    %             ur:           unit root case: I(0)|I(1)|MFI(1)|I(2). 
    %             n:            state dimension.
    %             s:            dimension of endogenous vars.
    %             y:            Txs matrix of observations of endogenous
    %                           vars.
    %             dt:           Txm matrix of obs for exogenous vars. 
    %             S:            integer; number of seasons per year. 
    %             time:         Tx1 vector of reals; time axis.
    %             
    properties
        % state space matrices
        theta = [];
        deviance = [];
        res = [];
        param = [];
        call = [];
        restrict = 0;
        Pbull = 0;
        urs = 0;
        ur = 'I(0)';
        n = 0;
        s = 0;
        y = [];
        dt = [];
        aic = Inf;
        bic = Inf;
        S = 1;
        time = [];
        V_theta = [];
        V_syst = [];
        betaquer = [];
    end
    
    methods
        % calculate log-likelihood (deviance)
        function result = cal_deviance(result)
            [qlike,tres] = cal_quasi_like(result.param,result.y,result.s,result.n,result.c,result.Pbull,result.det_res);
            result.deviance = qlike;
            result.res = tres; 
        end
        
            
        % generate impulse response function
        function ir = impulse(result,ML)
            
            if nargin<2
                ML = 100;
            end
            th = result.theta;
            if ~isa(th,'theta')
                error('impulse: can only be applied to a theta object!');
            end
            
            if strmatch(th.which,'poly')
                th = poly2ss_th(th);
            end
            p = size(th.C,1);
            m = th.m;
            
            ir = zeros(p,p+m,ML);
            ir(:,:,1)=eye(p);
            CAc = th.C;
            
            for j=1:(ML-1)
                ir(:,:,j+1)=CAc * [th.K,th.B];
                CAc = CAc * th.A;
            end
        end      

        function present(result,ylabs,which)
            s = result.s;
            n = result.n;

            % standard deviations.
            stds = sqrt(diag(result.V_syst));
            if nargin < 3
                which ='norm';
            end
            if nargin<2 
                for j=1:s
                    ylabs{j} = sprintf('y_%d',j);
                end
            end
            for j=1:n
                    xlabs{j} = sprintf('x_%d',j);
            end
            m = size(result.dt,2);
            for j=1:m
                    slabs{j} = sprintf('d_%d',j);
            end
            T = size(result.y,1);
            fprintf('%%%%%%%%%%%%%%%%%%%%%%%')
            fprintf(' Estimated unit root model. Unit root structure: %s',result.ur);
            fprintf('\n T: %d, s: %d, m: %d, n: %d, c: %d',T,result.s,m,result.n,result.urs)
            LL = -0.5*result.deviance;
            fprintf('\n Likelihood: %6.2f. Deviance: %6.2f',LL,-2*LL);
            npar = length(result.param);
            AIC = -2*LL+2*npar;
            BIC = -2*LL+log(T)*npar;
            fprintf('\n AIC:  %6.2f. BIC: %6.2f\n',AIC,BIC);
            fprintf('\n\n Omega (from final estimates)\n')
            disp(result.theta.Omega)
            fprintf('\n %%%%%%%%%%%%%%%%%%%%% ')
            zeros = sort(eig(result.theta.A-result.theta.K*result.theta.C));
            poles = sort(eig(result.theta.A)); 
            matprint(poles(:)',abs(poles(:)),{'Poles'});
            matprint(zeros(:)',abs(zeros(:)),{'Zeros'});
            fprintf('\n %%%%%%%%%%%%%%%%%%%%% ')

            fprintf('\n A: \n');
            matprint(result.theta.A,stds(1:n^2),xlabs,xlabs);
            fprintf('\n C: \n')
            matprint(result.theta.C,stds(n^2+[1:n*s]),ylabs,xlabs);
            fprintf('\n K: \n')
            matprint(result.theta.K,stds(n^2+n*s+[1:n*s]),xlabs,ylabs);
            fprintf('\n D: \n')
            matprint(result.theta.D,stds(n^2+2*n*s+[1:m*s]),ylabs,slabs);
       
            switch which
                case 'beta'
                    A = result.theta.A;
                    K = result.theta.K;
                    C = result.theta.C;
                    Abar = A-K*C;
                    ty = result.y - result.dt*result.theta.D';
                    r = n-result.urs;
                    if sum(strcmp(fieldnames(result), 'betaquer')) == 1
                        betaquer = result.betaquer;
                    else
                        betaquer = zeros(s,r);
                        betaquer(1:r,1:r)= eye(r);
                    end
                    [~,alphahat,betahat,~,V_alpha,V_beta] = RH_specm(ty,Abar,K,r,'y',5,betaquer);
                    % rlabs
                    for j=1:r
                        rlabs{j} = sprintf('F%d',j);
                    end
                    % alpha 
                    stdal = sqrt(diag(V_alpha));
                    fprintf('\n alpha: \n')
                    matprint(alphahat,stdal,ylabs,rlabs);

                    % beta
                    stdbeta = sqrt(diag(V_beta));
                    fprintf('\n beta: \n')
                    matprint(betahat,stdbeta,ylabs,rlabs);
                    fprintf('\n %%%%%%%%%%%%%%%%%%%%%\n ')
                otherwise
                    fprintf('\n %%%%%%%%%%%%%%%%%%%%%\n ')
            end
        end
    end
end
