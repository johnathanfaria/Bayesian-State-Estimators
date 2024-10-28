% Extended Kalman Filter / Unscented Kalman Filter / Moving Horizon
% Estimator applied to SONEGO et al.,2018
% EKF implementation according to Rawlings et al. (2019)
% UKF implementation according to Julier et al. (2000)/Wan et al. (2000)
% MHE implementation according to Rao & Rawlings (2002)
%
% by Gilson Campani and Johnathan Faria (2020)

clc;
clear all;
close all;

% odsets to derivation
opcoes = odeset('AbsTol',1e-6,'RelTol',1e-3);
% Numeric derivates conditions
dPar = 1e-8;
dx = 1e-8;
   

% Model parameters Kinetics
Yxs = 0.0415; s_Yxs = 0.0022;
Yes = 0.452; s_Yes = 0.006;
mi_max = 0.125; s_mimax = 0.002;
Ks = 25.1; s_Ks = 1.8;
Kis = 131.8; s_Kis = 9.3;
Ce_max = 86.1; s_Cemax = 1.7;   
n = 0.22; s_n = 0.03;
        
% Parameters array list
Par(1) = Yxs;
Par(2) = Yes;
Par(3) = mi_max;
Par(4) = Ks;
Par(5) = Kis; 
Par(6) = Ce_max;
Par(7) = n;


% Initial state conditions
% x = [Cx Cs Ce V]
% x = a posteriori // xminus = a priori // xk = real // xk_id = ideal

C0= [50 1e-8 1e-8 1.5];
% C0 = [40 1e-8 1e-8 1.8];
xk(:,1) = C0; % Real process with noise
% xk_id(:,1) = xk(:,1);% Ideal process without noise

% Sample time
t0=0;
tf=13;
tAmos=1e-;
tspamAmos=[t0:tAmos:tf];
imax = length(tspamAmos)-1;
tAmos_id=0.1;
tspamAmos_id=[t0:tAmos_id:tf];
Dt = tAmos;
Dt_id = tAmos_id;


%Plots COV
xk_up = xk;
xk_down = xk;

% Measurament parameters
C = [1 0 0 0;
    0 1 0 0;
    0 0 0 1];

a = 0.2;
b = 0.1;
c = 0.2;

cM_aux = [ a*xk(1);
           b*xk(2)
           c*xk(4)]; 
       
R = diag(cM_aux.^2); 
yk = C*xk + mvnrnd(zeros((length(R)),1),R)';

% Covariance matrix of paramters
Q = diag([s_Yxs,s_Yes,s_mimax,s_Ks,s_Kis,s_Cemax,s_n].^2);
% sqrt_Q = [s_Yxs,s_Yes,s_mimax,s_Ks,s_Kis,s_Cemax,s_n,s_ke,s_kw];

% Initial state and covariance matrix and measurement
Pminus = eye(4);
% xminus(:,1) =  C0;
xminus(:,1) =  [40 1e-8 1e-8 1.80];   

% A-posteriori state and covariance matrix
x(:,1) = xminus(:,1);
P = Pminus;

% MHE parameters e initial variables
HN = 1; % hozizon legth 6 - it has to be >=2 to avoid problem in initial guess X0 of function MHE()
alfa = 0; % forgetting factor
Pminus0 = Pminus;
Pminus_aux = Pminus0; % Pminus(T-N) used by function MHE

    
% UKF parameters e initial variables
xa = [x(:,1);zeros(length(diag(Q))+length(diag(R)),1)];% augmented vector: xa = [x;w;v]
Pa = blkdiag(P,Q,R); % augmented covariance matrix: xa = [x;w;v]
n = length(x(:,1)); % state dimension
na = length(xa(:,1)); % augmented state dimension
nS = 2*na+1;% number of sigma points
Sigma = zeros(nS,na); % initialization of matrix with sigma points
SigmaMinus = zeros(nS,na); % initialization of propagated matrix with sigma points

% definition of parameter kapa: 
% 0 < kapa + na < 3 = kurtosis of Gaussian distribution 
% (Romanenko et al. 2004)

% Heuristics from Vachhani et. al. (2006) 
% This Value Guarantees that Pa is Positive Definite
kapa = 1;

% % Heuristics from Kandepu et. al. (2008)
% alpha = 0.1;      % alpha Should be a Small Number Between 0 and 1
% kapa = alpha^2*na - na;

% % Heuristics from Julier et al. (2000): 0 < kapa + 3 < 3 = kurtosis of Gaussian distribution
% kapa = 3 - na; 

% Trial and Error
% kapa = -6;         % Trial and Error - Lima

W(1) = kapa/(na+kapa); % initialization of weight for sigma points
W(2:nS) = 1/(2*(na+kapa)); % initialization of weight for sigma points

tspam_aux = linspace(t0,tAmos,5);% Subdivision of tspam

% Process simulation without noise
% [T,X] = ode45(@edo1,tspam_aux,xk_id(:,1),opcoes,Par,F,ke,kw);% Process simulation
% xk_id = X';       

% e_plot(1,:) = 0;
% Par_aux = Par + mvnrnd(zeros(1,9),Q);
% Par_aux = Par + mvnrnd(zeros((length(Q)),1),Q);
% GQ_aux = zeros(9,1);

tstart = tic;
% State estimator method ('EKF', 'UKF' or 'MHE')
F = 0.56;
method = 'MHE';
PDT_i(:,1) = 0;

for i=1:imax
     
   %%%%%%% comeco do stripping
    if xk(3,end) <= 34.18
        ke = 0;
        kw = 0;
    else
        ke = 0.0656; s_ke = 0.0016;
        kw = 0.00443; s_kw = 0.00002;
        PDT_i(:,end+1) = xk(3,end)*xk(4,end);
       
    end

    % Real process simulation
%     Par_aux = Par;
    Par_aux = Par + mvnrnd(zeros((length(Q)),1),Q);
    tspam_aux = linspace(tspamAmos(i),tspamAmos(i)+Dt,5);% Subdivision of tspam
    [T,X] = ode45(@edo1,tspam_aux,xk(:,i),opcoes,Par_aux,F,ke,kw);% Process simulation with plant mismatch
    xk(:,i+1) = X(end,:)';
    xk_aux = X(end,:)';
%     tspam_aux(1)
    
    % Measurement
    yk(:,i+1) = C*xk(:,i+1) + mvnrnd(zeros((length(R)),1),R)'; % Normal noise addition to measurement
%     yk(:,i+1) = C*xk(:,i+1);
    cM_aux = [a*xk(1,end);
                b*xk(3,end)
                c*xk(4,end)];
    R = diag(cM_aux.^2);
    
    if rem(i,500) == 0 
        jjj = i/500;
        ty_plot(jjj+1) = (i/1000);
        y_plot(:,jjj+1) = yk(:,end);
        
    end
    
    
    % Estimation base on option (EKF, UKF or MHE)
    switch method
        case 'EKF'
            % State prediction
            xk_aux = x(:,i);
            tspam_aux = linspace(tspamAmos(i),tspamAmos(i)+Dt,5);% Subdivision of tspam
            [T,X] = ode45(@edo1,tspam_aux,xk_aux,opcoes,Par,F,ke,kw);% Process simulation
 
            % Jacobinan matrix A <-> x =~ A*x + G*Q
            for j=1:length(xk_aux)
                xk_aux_inc = xk_aux;
                xk_aux_inc(j) = xk_aux(j)*(1+dx);
                [T,Xinc] = ode45(@edo1,tspam_aux,xk_aux_inc,opcoes,Par,F,ke,kw);% Process simulation with x increment
                A(:,j) = [(Xinc(end,1)-X(end,1))/(xk_aux_inc(j));(Xinc(end,2)-X(end,2))/(xk_aux_inc(j));(Xinc(end,3)-X(end,3))/(xk_aux_inc(j));(Xinc(end,4)-X(end,4))/(xk_aux_inc(j))];
            end
            
            % Jacobinan matrix G <-> x =~ A*x + G*Q
            for j=1:length(Par)
                Par_inc = Par;
                Par_inc(j) = Par(j)*(1+dPar);
                [T,Xinc] = ode45(@edo1,tspam_aux,xk_aux,opcoes,Par_inc,F,ke,kw);% Process simulation with Par increment
                G(:,j) = [(Xinc(end,1)-X(end,1))/(Par_inc(j));(Xinc(end,2)-X(end,2))/(Par_inc(j));(Xinc(end,3)-X(end,3))/(Par_inc(j));(Xinc(end,4)-X(end,4))/(Par_inc(j))];
            end
            
            % State and covariance forecasting
            [T,Xminus] = ode45(@edo1,tspam_aux,x(:,i),opcoes,Par,F,ke,kw);
            xminus(:,i+1) = Xminus(end,:)';% Simulated state
            Pminus = A*P*A' + G*Q*G';   
            
            % A posteriori state and covariance ("correction")
            L = Pminus*C'/(C*Pminus*C'+R);
            x(:,i+1) = xminus(:,i+1) + L*(yk(:,i+1)-C*xminus(:,i+1));
            P = Pminus - L*C*Pminus;

            
        case 'UKF'

            % First sigma point calculation
            sqrtP = real(sqrtm((na+kapa)*Pa)); % real() to avoid imaginary
            Sigma(1,:) = xa;
            
            % First sigma point propagation
            [T,X] = ode45(@edo1,tspam_aux,Sigma(1,1:n),opcoes,Par_aux,F,ke,kw);
            w = Sigma(1,n+1:n+length(diag(Q))); % process noise from sigma points
            v = Sigma(1,n+length(diag(Q))+1:end); % measurement noise from sigma points
            Xaux = X(end,:);% state prediciton with sigma's process noise
            SigmaMinus(1,:) = [Xaux  w  v];
            
            for j=1:na
                % Further sigma points calculation
                Sigma(j+1,:) = xa' + sqrtP(j,:);
                Sigma(j+1+na,:) = xa' - sqrtP(j,:);
                
                % Sigma points time propagation (2:na)
                [T,X] = ode45(@edo1,tspam_aux,Sigma(j+1,1:n),opcoes,Par_aux,F,ke,kw);
                w = Sigma(j+1,n+1:n+length(diag(Q))); % process noise from sigma points
                v = Sigma(j+1,n+length(diag(Q))+1:end); % measurement noise from sigma points
                Xaux = X(end,:);% state prediciton with sigma's process noise
                SigmaMinus(j+1,:) = [Xaux  w  v];
                % Sigma points time propagation (na+1:2na+1)                
                [T,X] = ode45(@edo1,tspam_aux,Sigma(j+1+na,1:n),opcoes,Par_aux,F,ke,kw);%(tspam_aux,Sigma(j+1,:));
                w = Sigma(j+1+na,n+1:n+length(diag(Q))); % process noise from sigma points
                v = Sigma(j+1+na,n+length(diag(Q))+1:end); % measurement noise from sigma points
                Xaux = X(end,:);% state prediciton with sigma's process noise
                SigmaMinus(j+1+na,:) = [Xaux  w  v];
            end
            
            % Sigma measurements calculation
            SigmaY = (C*SigmaMinus(:,1:n)')' + Sigma(:,n+length(diag(Q))+1:end); % noise addition
            yminus = (W*SigmaY)';
            
            % A priori mean and covariance calculation
            xminus(:,i+1) = (W*SigmaMinus(:,1: n))';
            Pminus = (SigmaMinus(:,1:n)' - xminus(:,i+1))*diag(W)*(SigmaMinus(:,1:n)' - xminus(:,i+1))';

            % A posteriori calculation
            Pyy = (SigmaY' - yminus)*diag(W)*(SigmaY' - yminus)';
            Pxy = (SigmaMinus(:,1:n)' - xminus(:,i+1))*diag(W)*(SigmaY' - yminus)';
            L = Pxy/Pyy; % Kalman gain
            x(:,i+1) = xminus(:,i+1) + L*(yk(:,i+1)-C*xminus(:,i+1)); % a posteriori state estimate
            xa = [x(:,i+1);zeros(length(diag(Q))+length(diag(R)),1)];
            P = Pminus - L*Pyy*L'; % a posteriori covariance estimate
            Pa = blkdiag(P,Q,R);
                   
        case 'MHE'
            % Jacobinan matrix G <-> x =~ A*x + G*Q
            xk_aux = x(:,i);
            for j=1:length(Par)
                Par_inc = Par;
                Par_inc(j) = Par(j)*(1+dPar);
                [T,Xinc] = ode45(@edo1,tspam_aux,xk_aux,opcoes,Par_inc,F,ke,kw);% Process simulation with Par increment
                G(:,j) = [(Xinc(end,1)-X (end,1))/(Par_inc(j));(Xinc(end,2)-X(end,2))/(Par_inc(j));(Xinc(end,3)-X(end,3))/(Par_inc(j));(Xinc(end,4)-X(end,4))/(Par_inc(j))];
            end
            alfa1 = 1e-14;
            invQ = inv(G*Q*G'+ alfa1.*eye(4));
            invR = inv(R); 
           % MHE calculation
            [xopt,Pminus_aux] = MHE(tspamAmos,i,Pminus0,Pminus_aux,invR,invQ,C,HN,alfa,yk,x,R,G,Q,Par,F,ke,kw);
            x(:,i+1) = xopt; 
            
        case 'COV'
            
           [T,X] = ode45(@edo1,tspam_aux,xk(:,i),opcoes,Par,F,ke,kw);
           xk(:,i+1) = X(end,:)';
            xk_up(:,i+1) = X(end,:)';
            xk_down(:,i+1) = X(end,:)';
            mi_t(:,i+1) = MI(xk(:,i),Par);
            for j=1:length(Par)
                Par_inc = Par;
                Par_inc(j) = Par(j)*(1+dPar);
                [T,Xinc] = ode45(@edo1,tspam_aux,xk(:,i+1),opcoes,Par_inc,F,ke,kw);% Process simulation with Par increment
                x_cov(:,j) = [(Xinc(end,1)-X(end,1))/(Par_inc(j));(Xinc(end,2)-X(end,2))/(Par_inc(j));(Xinc(end,3)-X(end,3))/(Par_inc(j));(Xinc(end,4)-X(end,4))/(Par_inc(j))];
                
            end
            error = 2.034515*sqrt(diag(x_cov*Q*x_cov'));    
            xk_up(:,i+1) = xk_up(:,i+1) + error;
            xk_down(:,i+1) = xk_down(:,i+1) - error;
           
    end
        
         
     
     e = [(x(1,end)-xk(1,end)),(x(2,end)-xk(2,end)),(x(3,end)-xk(3,end)),(x(4,end)-xk(4,end))];
     e_plot(i+1,:) = (e');
     e2(:,i+1) = sum(e.^2/xk(:,end)');
        
        
% % % %         xk_up(:,i) = xk_up(:,i) + error;
% % % %         xk_down(:,i) = xk_down(:,i) - error;
% % % %         

    if i >=6250
         F = 0;
    end
end
% % % 
% % % temp = [tspamAmos'; xk(1,:)'; x(1,:)';xk(2,:)'; x(2,:)'; xk(3,:)'; x(3,:)'];
% % % concentracoes0=strcat('concentracoes', method);
% % % concentracoes=[concentracoes0 '.xls'];
% % % save (concentracoes,    'temp',  '-ASCII')
% % % 
% % % temp1 = [ty_plot';y_plot(1,:)';y_plot(2,:)';y_plot(3,:)'];
% % % medidas=strcat('concentracoes_medidas', method);
% % % medidas=[medidas '.xls'];
% % % save (medidas,    'temp1',  '-ASCII')


disp(method);
toc(tstart);



jmax = length(PDT_i)-1;
soma = 0;
for j=2:jmax
    if j == 2 || j == jmax
        soma = soma + PDT_i(j);
    elseif rem(j,2) == 1
        soma = soma + (2*PDT_i(j));
    elseif rem(j,2) == 0
        soma = soma + (4*PDT_i(j));
    end
end


PDT_I = (Dt/3)*soma*ke;
CET = ((xk(3,end)*xk(4,end))+ PDT_I);
disp("Concentração total de etanol");
disp(CET);


% e2 = [((x(1,:)-xk(1,:))),((x(2,:)-xk(2,:))),((x(3,:)-xk(3,:))),((x(4,:)-xk(4,:)))]./[xk(1,:),xk(2,:),xk(3,:),xk(4,:)];
ex1 = 100*((x(1,:)-xk(1,:)))./xk(1,:);
ex2 = 100*((x(2,:)-xk(2,:)))./xk(2,:);
ex3 = 100*((x(3,:)-xk(3,:)))./xk(3,:);
ex2(1)=ex2(2);
ex3(1)=ex3(2);

disp('Residual Sum of Squares per tspam');
EMQ = sqrt(sum(e2)/length(e2));
disp(EMQ);
       


figure(1)
plot(tspamAmos,xk(1,:)','-k',tspamAmos,xk(2,:),'-b',tspamAmos,xk(3,:)','-r',...
    tspamAmos,x(1,:)','-.k', tspamAmos,x(2,:),'-.b',tspamAmos,x(3,:)','-.r',...
    ty_plot,y_plot(1,:),'.k',ty_plot,y_plot(2,:),'.b');
   
xlabel('Time [h]');
ylabel('Concentration [gL^{-1}]')
title('a)UKF')

ylim([0 140]);
legend('C_X','C_S','C_E',...
    'C_X estimated','C_S estimated','C_E estimated',...
    'C_X measured','C_S measured')
% figure (2)
% plot(tspamAmos,e_plot(:,1),'-k',tspamAmos,e_plot(:,2),'-b',tspamAmos,e_plot(:,3),'-r')%,tspamAmos,e_plot(:,4),'-g')
% xlabel('Time [h]')
% ylabel('Error')
% % title('E')
% ylim([-10 10]);
% legend('C_X','C_S','C_E');

figure(2)
plot(tspamAmos,ex1,'-k',tspamAmos,ex2,'-b',tspamAmos,ex3,'-r');
xlabel('Time [h]','FontSize', 14)
ylabel('Error [%]','FontSize', 14)
ylim([-80 40]);
legend('C_X','C_S','C_E');


% figure(4)
% plot(tspamAmos,x(2,:))
% 
% figure(3)
% plot(tspamAmos,xk(1,:)','-k',tspamAmos,xk(2,:)','-b',tspamAmos,xk(3,:)','-r',tspamAmos,xk(4,:)','-g',...
%     tspamAmos,xk_up(1,:)','-.k',tspamAmos,xk_up(2,:)','-.b',tspamAmos,xk_up(3,:)','-.r',tspamAmos,xk_up(4,:)','-.g',...
%     tspamAmos,xk_down(1,:)','-.k',tspamAmos,xk_down(2,:)','-.b',tspamAmos,xk_down(3,:)','-.r',tspamAmos,xk_down(4,:)','-.g')
% xlabel('Tempo [h]')
% ylabel('Concentração [gL^{-1}]')
% legend('C_X - Modelo','C_S - Modelo ', ' C_E - Modelo','V - Modelo',...
%     'C_X - Intervalo de Confiança','C_S - Intervalo de Confiança ', ' C_E - Intervalo de Confiança','V - Intervalo de Confiança')
% 
% yyaxis right
% % plot(tspamAmos,xk(4,:)','-g',tspamAmos,xk_up(4,:)','-.g',tspamAmos,xk_down(4,:)','-.g')
% ylim([0 140])
% ylabel('Volume [L]')
% legend('C_X - Modelo','C_S - Modelo ', ' C_E - Modelo','V - Modelo',...
%     'C_X - Intervalo de Confiança','C_S - Intervalo de Confiança ', ' C_E - Intervalo de Confiança','V - Intervalo de Confiança')
%     


   