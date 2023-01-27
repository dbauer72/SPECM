classdef theta_urs
    % Ljung's theta object to save ARMAX or state space systems.
    % Fields: which: 'poly': (a,b,d). 
    %         which: 'SS': (A,B,C,D,K)
    %         Omega ... sxs innovation covariance matrix.
    %         nd    ... integer;
    %         m     ... integer; dimension of exogenous inputs. 
    % special object for unit root systems
    properties
        % state space matrices
        A = [];
        B = [];
        C = [];
        D = [];
        K = [];
        % polynomial representation
        a = 1;
        b = 1;
        d = [];
        % innovation noise variance
        Omega = 1;
        % flag pointing to active representation.
        which {mustBeMember(which,{'SS','poly'})} = 'SS'
        % 
        nd = 0;
        m = 0;
        ur = '';
        urs = [];
        num_param = 0;
    end
    
    methods
        % check dimensions 
        function ok = check_dims(th)
            ok = 1
            if strmatch(th.which,'poly') % polynomials are specified
                [ra,ca] = size(th.a);
                [rb,cb] = size(th.b);
                [rd,cd] = size(th.d);
                if (ra ~= rb) 
                    ok = 0;
                    disp('Dimensions of a and b do not match!')
                end
                if (ra ~= rd) 
                    ok = 0;
                    disp('Dimensions of a and d do not match!')
                end
                if (rem(ca,ra) ~= 0) 
                    ok = 0;
                    disp('Dimensions of a do not match!')
                end
                if (rem(cb,ra) ~= 0) 
                    ok = 0;
                    disp('Dimensions of b do not match!')
                end                
                if (rem(cd,th.m) ~= 0) 
                    ok = 0;
                    disp('Dimensions of d do not match!')
                    th.d = zeros(ra,th.m);
                end

            else % it is a state space system
                [ra,ca] = size(th.A);
                [rb,cb] = size(th.B);
                [rc,cc] = size(th.C);
                [rd,cd] = size(th.D);
                [rk,ck] = size(th.K);
                if (ra ~= rb) || (ra ~= rk) 
                    ok = 0;
                    disp('System order not uniform.')
                end
                if (ra ~= ca)
                    ok = 0;
                    disp('A is not square!')
                end
                if (rc ~= rd)
                    ok = 0
                    disp('C and D do not have the same number of rows!')
                end
                if (cd ~= cb)
                    ok = 0;
                    disp('B and D do not have the same number of columns!')
                end
                if (rc~= ck)
                    ok = 0;
                    disp('The number of rows of C does not match the number of columns of K!')
                end


            end
            [ro,co] = size(th.Omega);
            if (ro ~= co)
                ok = 0;
                disp('Dimension of Omega does not match!')
            end
        end



        % generate impulse response function
        function ir = impulse(th,ML);
            
            if nargin<2
                ML = 100;
            end;
            if ~isa(th,'theta_urs')
                error('impulse: can only be applied to a theta object!');
            end;
            
            if strmatch(th.which,'poly')
                th = poly2ss_th(th);
            end;
            p = size(th.C,1);
            m = size(th.D,2);
            
            ir = zeros(p,p+m,ML);
            ir(:,:,1)=[eye(p),th.D];
            CAc = th.C;
            
            for j=1:(ML-1)
                ir(:,:,j+1)=CAc * [th.K,th.B];
                CAc = CAc * th.A;
            end;
        end;
        
        % generate bode plot
        function [mod_tf,phase_tf,tf] = bode(th,fr);
            if nargin<2
                fr = exp(([1:512]/51.2-10)+ log(pi));
            end;
            
            switch th.which
                case 'SS'
                    [p,n] = size(th.C);
                    m = th.m;
                    tf = zeros(p,p,length(fr));
                    for j=1:length(fr)
                        z = exp(sqrt(-1)*fr(j)); % exp(i omega)
                        tf(:,:,j) = eye(p) + z*th.C*inv(eye(n) - z*th.A)*th.K;
                    end;
                    mod_tf = abs(tf);
                    phase_tf = angle(tf);
                case 'poly'
                    p = size(th.a,1);
                    m = th.m;
                    tf = zeros(p,p,length(fr));
                    a = th.a;
                    na = length(a);
                    b = th.b;
                    nb = length(b);
                    if na>nb
                        b = [b,zeros(1,na-nb)];
                    else
                        a = [a,zeros(1,nb-na)];
                    end;
                    for j=1:length(fr)
                        z = exp(sqrt(-1)*fr(j));
                        zz = z.^[0:max(na,nb)-1];
                        tf(:,:,j) = (zz*a(:))\(zz*b(:));
                    end;
                    mod_tf = abs(tf);
                    phase_tf = angle(tf);
            end;
        end;
        
        % return polynomials
        function [th,a,b,d] = th2poly(th);
            th = ss2poly(th);
            
            a = th.a;
            b = th.b;
            d = th.d;
        end;
        
        % return state space representation
        function [th,A,B,C,D,K] = th2ss(th);
            th = poly2ss_th(th); 
            A = th.A;
            B= th.B;
            C = th.C;
            D = th.D;
            K = th.K;           
        end;
        
        % difference transfer function. 
        function th = diff_tf(th);
            [th,A,B,C,D,K]=th2ss(th);
            [n,s] = size(K);

            Ad = [zeros(s,s),-C;zeros(n,s),A];
            Kd = [-eye(s);K];
            Cd = [eye(s),C];

            n=n+s;
            th = ss2th(Ad,B,Cd,D,Kd,th.Omega);
        end;
        
        % take inverse 
        function thi = inverse_TF(th);
            % calculates the inverse TF, but only of k!
            [th,A,B,C,D,K]=th2ss(th);
            
            thi = th;
            thi.A = th.A-th.K*th.C;
            thi.K= th.K;
            thi.C = -th.C;
        end
    end
end
