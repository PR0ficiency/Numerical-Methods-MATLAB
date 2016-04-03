% Craig Medlin
% Feb. 15, 2016
% AERO 220
% HW 2

% Use a 'switch' statement to allow selection of different problems

%problem = 6 ; % choose which problem you want to run

problem = input('Choose a problem (1-6): '); % Or ask the user


switch(problem)
    
    case 1
        %% Problem 1: Find steady-state concentration via fixed-point
        fprintf('----------- Problem 1 -----------\n');
        clear all;
        % Initial Function  ->  V(dc/dt) = W - Qc-kV*SQRT(c)
        
        g1 = @(c) ((10-c)/2.5)^2;
        g2 = @(c) (10-2.5*(c^0.5));
        
        % Given values
        V = 10^6; % m^3
        Q = 10^5; % m^3/yr
        W = 10^6; % g/yr
        k = 0.25; % g^0.5/(m^1.5 * yr)
        
        eps=0.1;
        
        for i=1:2
            %Initial Guess
            c0 = 4; % g/m^3
            
            switch i
                case 1
                    g = g1;
                case 2
                    g = g2;
                otherwise
                    fprintf('ERROR: Bounds on function');
            end
            
            while abs(g(c0)-c0)>eps
                c0 = g(c0);
            end
            
            g
            fprintf('Converges for g%i at %f\n',i,c0);
            
            
        end
        
        
        
    case 2
        %% Problem 2
        fprintf('----------- Problem 2 -----------\n');
        f = @(x) -x^3 + 2*x - exp(-x) + 3;
        x0 = 2.3;
        N = 5;
        
        syms x;
        dfx = diff(f(x));
        
        % Part A
        f_eval = f(x0);
        dfx_eval = vpa(subs(dfx,x,x0));
        
        X = [x0];
        ANS = [];
        
        if f_eval ~= 0 && dfx_eval ~= 0
            for i=1:5
                x1 = x0;
                x0 = x0 - f(x0)/vpa(subs(dfx,x,x0));
                
                X = [X; double(x0)];
            end
            
            for i=1:5
                en = abs(X(i)-X(6));
                lnen = log(en);
                
                ANS = [ANS; i-1 X(i) en lnen];
            
            end
            
            ANS = [ANS; 5 X(6) abs(X(6)-X(6)) log(0)];
            
            alpha = (ANS(5,4)-ANS(4,4))/(ANS(4,4)-ANS(3,4));
            lambda = ANS(4,3)/(ANS(3,3)^alpha);
            
            names = {'Iteration', 'Estimated_Root', 'Absolute_Error', 'lnen'};
            array2table(ANS, 'VariableNames', names)
            
            fprintf('\t\t  Estimated Root: %f\n\tOrder of Convergence: %f\nAsymptotic Error Constant: %f\n\n',ANS(6,2),alpha,lambda);
            
        end
        
        % Part B
        a = -3;
        b = 3;
        
        flag = 0;
        Xi = [];
        Yi = [];
        n = 1;
        for i=a:0.1:b
            Xi = [Xi;n];
            Yi = [Yi;abs(vpa(subs(dfx,x,i)))];
            if  Yi(n)>= 1
                flag = 1;
            end
            n = n+1;
        end
        
        if flag == 0
            fprintf('Analytic evaluation shows Newton''s method converges over [%i,%i].\n', a, b);
        else
            fprintf('Analytic evaluation shows Netwon''s method does NOT converge over [%i,%i].\n', a, b);
        end
        
        ezplot(abs(dfx),[-3 3]);
        title('|g''(x)| for x \in [-3, 3]');
        xlabel('X');
        ylabel('Y');
        
    case 4
        %% Problem 4
        fprintf('----------- Problem 4 -----------\n');
        A = [2 1 -1; 4 1 2; 6 1 1];
        b = [-2; 4; 6];        
        
        [x, A2,B2] = GaussElim(A, b);
        
        array2table(x,'VariableNames', {'x1', 'x2', 'x3'})
        
    case 5
        %% Problem 5 - Operate on matricies
        fprintf('----------- Problem 5 -----------\n');
        A = [-3 1 -2; 2 3 0; -1 2 3];
        B = [3 -2 5; 2 -4 1; -4 1 6];
        
        % Part A
        BA = B*A, B_Cubed = B^3,AA_T = A*A'
        
        % Part B
        detA = mDeter2(A)
        detB = mDeter2(B)
        
        % Part C
        [L_A U_A] = LUdec(A)
        [L_B U_B] = LUdec(B)
        
        % Part D
        
        
        
    otherwise
        error('Invalid Problem Number');
end
    
