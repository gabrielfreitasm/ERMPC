clear, clc, close all
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

% classical angular positioning system adapted from Kwakernaak and Sivan
% (1972) 
% Matrizes politopicas
A1 = [1 0.1; 0 0.99];
A2 = [1 0.1;0 0];
B = [0;0.1*0.787];

% Project restrictions
umax=2;

% Parameters Qc Rc
Qc = [1 0; 0 0]; Rc = 2e-5;

% ONLINE Model
% YALMIP/OPTIMIZER/SEDUMI Optimization
n=size(B,1); m=size(B,2);
Q=sdpvar(n,n,'symmetric');
X=sdpvar(m,m,'symmetric');
Y=sdpvar(m,n,'full');
gama=sdpvar(1);

xk = sdpvar(2,1);
% M1 Matrix - Restriction considering vertix 1
m11=[Q (A1*Q+B*Y)' (sqrt(Qc)*Q)' (sqrt(Rc)*Y)'];
m12=[(A1*Q+B*Y) Q zeros(n) zeros(n,m)];
m3=[(sqrt(Qc)*Q) zeros(n) gama*eye(n) zeros(n,m)];
m4=[(sqrt(Rc)*Y) zeros(n,m)' zeros(n,m)' gama*eye(m)];
M1=[m11;m12;m3;m4]; % M1

% M2 Matrix - Restriction considering vertix 2
m21=[Q (A2*Q+B*Y)' (sqrt(Qc)*Q)' (sqrt(Rc)*Y)'];
m22=[(A2*Q+B*Y) Q zeros(n) zeros(n,m)];
M2=[m21;m22;m3;m4]; % M2

LMIs=[Q >= 0, M1 >= 0, M2>=0, X>=0, gama >= 0]; % LMIs restrictions
LMIs=[LMIs, [1 xk';xk Q]>=0]; % assures that J is bounded
LMIs=[LMIs, [X Y;Y' Q]>=0, X<=umax^2]; % adds a control signal restriction

ops = sdpsettings('solver','sedumi','sedumi.eps',1e-5); % sets the solver and its accuracy
model = optimizer(LMIs, gama,ops,xk,{Y,Q});

% Initial conditions
x = [1;0];
xs = [1;0];
QY = model{x}; % returns Y and Q for the initial condition x
Fs = QY{1}/(QY{2}); % static model

tsim = 20; ts=0.1;
num = ceil(tsim/ts);
alpha = zeros(1,num);
u = alpha;
us = alpha;
F = cell(1,num); % dynamic gain - calculated in each iteraction

for k = 1:num

    alpha(k)=10*rand; % random number in (1,10)
    A = [1 0.1; 0 1-0.1*alpha(k)]; % updates the time varying model
    QY = model{x(:,k)}; % updates the Y and Q decision variables
    F{k} = QY{1}/(QY{2}); % is equal to the static gain for k = 1

    u(k) = F{k}*(x(:,k));
    x(:,k+1) = A*x(:,k) + B*u(k);
    us(k) = Fs*(xs(:,k));
    xs(:,k+1) = A*xs(:,k) + B*us(k);

end

t = 0:ts:tsim;

% figures
figure(1)
subplot(1,2,1)
plot(t,x(1,:),'k', t,xs(1,:),'b--','LineWidth',2); % plots only the first state, since theta = x(1,:)
xlabel('Time (sec)','FontSize',13)
ylabel('$\theta$ (rad)','FontSize',13)

title("Angular position")
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',13)
grid on
subplot(1,2,2)
plot(t(1:end-1),u,'k-',t(1:end-1),us,'b--','LineWidth',2);
xlabel('Time (sec)','FontSize',13)
ylabel('u (volts)','FontSize',13)
title("Control Signal")
grid on
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',13)
legend('ERMPC','Static MPC','FontSize',13)





