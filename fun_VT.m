
% Cost function for MHE optimization
% 
% Gilson - set/2017 and Johnathan dez/2020

function [VT] = fun_VT(x,dt,invQ,invR,C,alfa,invP,y,xT_N,Par,F,ke,kw)
opcoes = odeset('AbsTol',1e-6,'RelTol',1e-3);

% System mathematical model
% x = [Cx Cs Ce V]

% VTw total cost through horizon (model deviation)
for i=1:(length(y(1,:))-1)
    % w vector calculation
    [T,Ax] = ode45(@edo1,[0;dt],x(:,i),opcoes,Par,F,ke,kw);
    Ax = Ax(end,:)';
    W(i) = (x(:,i+1)-Ax)'*invQ*(x(:,i+1)-Ax);
end
VTw = sum(W);

% VTu total cost through horizon (measurement deviation)
for i=1:(length(y(1,:)))
     U(i) = (y(:,i)-C*x(:,i))'*invR*(y(:,i)-C*x(:,i));
end
VTu = sum(U);

% Arrival cost calculation
VTa = alfa*(x(:,1)-xT_N)'*invP*(x(:,1)-xT_N);

% Total cost

VT = VTw + VTu + VTa;
VT = VT;


end

