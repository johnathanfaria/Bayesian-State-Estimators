function [xopt,Pminus_aux] = MHE(tspamAmos,i,Pminus0,Pminus,invR,invQ,C,HN,alfa,yk,x,R,G,Q,Par,F,ke,kw)

% Moving Horizon Estimation - function
%
% Implementation according to Rao & Rawlings 2002
%
% Aug - 2017 Gilson Campani and edited in Dez - 2020  by Johnathan Faria

opcoes = odeset('AbsTol',1e-6,'RelTol',1e-3);

% options = optimoptions('fmincon');
options.MaxIterations = 8000;
options.OptimalityTolerance = 1.0000e-08;
% options.Algorithm = 'sqp';
options.StepTolerance = 1e-8;
options.OptimalityTolerance = 1e-10;
options.FunctionTolerance = 1e-10;
options.MaxFunctionEvaluations = 8000;
% options.Display = 'iter';
% options.MaxFunEvals = 6000;
% options = optimset('MaxFunEvals',3000,'MaxIter',3000);
options.Display = 'none';

% Sample time
tAmos = tspamAmos(2)-tspamAmos(1);
% Paramaters for numerical derivatives
dx = [1e-8;1e-8;1e-8;1e-8];

% MHE calculation
if i<=HN
    % Initial guess for optimization
    % State forecasting
    [T,X] = ode45(@edo1,tspamAmos(i:i+1),x(:,i),opcoes,Par,F,ke,kw);% Process simulation
    X0 = [x(:,[1:i]) X(end,:)'];
%     X0 = xk_id(:,[1:i+1]).*(1+0.001*randn(size(xk_id(:,[1:i+1]))));

    % Optimization
    xT_N = x(:,1);
    fun = @(x)fun_VT(x,tAmos,invQ,invR,C,alfa,inv(Pminus0),yk(:,[1:i+1]),xT_N,Par,F,ke,kw);
%     xopt = fmincon(fun,X0,[],[],[],[],zeros(size(X0)),[],[],options);
    xopt = fminunc(fun,X0,options);
%     xopt = fminsearch(fun,X0,options);
    xopt = xopt(:,end);
    Pminus_aux = Pminus0;% calculated Pminus
elseif i>HN
    % Initial guess for optimization
    % State forecasting
    [T,X] = ode45(@edo1,tspamAmos(i:i+1),x(:,i),opcoes,Par,F,ke,kw);% Process simulation
    X0 = [x(:,[i-HN+1:i]) X(end,:)'];
%     X0 = xk_id(:,[i-HN+1:i+1]).*(1+0.001*randn(size(xk_id(:,[i-HN+1:i+1]))));
    
    % Determination of x(T-N)
%     [T,X] = simulacaoX(tspamAmos(i-HN:i-HN+1),x(:,i-HN));% Process simulation
    xT_N = x(:,i-HN+1);
    
    % State simulation whithout w or r increment
    [T,X] = ode45(@edo1,tspamAmos(i-HN+1:i-HN+2),x(:,i-HN+1),opcoes,Par,F,ke,kw);% Process simulation
    % Jacobinan matrix A <-> x =~ A*x + G*w
    for j=1:length(x(:,i-HN+1))
        x_aux = x(:,i-HN+1);
        x_inc = x_aux;
        x_inc(j) = x_aux(j)+dx(j);
        [T,Xinc] = ode45(@edo1,tspamAmos(i-HN+1:i-HN+2),x_inc,opcoes,Par,F,ke,kw);% Process simulation with x increment
        A(:,j) = [(Xinc(end,1)-X(end,1));(Xinc(end,2)-X(end,2));(Xinc(end,3)-X(end,3));(Xinc(end,4)-X(end,4))]/dx(j);
    end
    
    % Covariance forecasting and MHE optimization
    Pminus = G*Q*G' + A*(Pminus-(Pminus*C'/(R+C*Pminus*C'))*C*Pminus)*A';%C*R*C'
    L = Pminus*C'/(C*Pminus*C'+R);%C*R*C'
    P = Pminus - L*C*Pminus;
    fun = @(x)fun_VT(x,tAmos,invQ,invR,C,alfa,inv(P),yk(:,[i-HN+1:i+1]),xT_N,Par,F,ke,kw); %yk(:,[i-HN+1:i+1])
%     xopt = fmincon(fun,X0,[],[],[],[],zeros(size(X0)),[],[],options);
    xopt = fminunc(fun,X0,options);
%     xopt = fminsearch(fun,X0,options);
    xopt = xopt(:,end);         
    Pminus_aux = Pminus;% calculated Pminus
end
    
end
