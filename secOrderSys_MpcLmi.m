clear, clc, close all
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

% classical angular positioning system adapted from Kwakernaak and Sivan
% (1972) 
% Polytopic matrices
%C:\Users\aaron\OneDrive\Área de Trabalho\Backup do Drive\2022.1\LMI\Programas\Class 5 Polytopic Uncertainties\code5_continuous_polytopic.m
A1 = [-4 4;-5 0] + [-2 2;-1 4];
A2 = [-4 4;-5 0] - [-2 2;-1 4];
B = [1;2];
C = [1 1];
D = 0;

% Project restrictions
umax=4;

% Parameters Qc Rc
Qc = [1 0; 0 0]; Rc = 2e-5;

% ONLINE Model
% YALMIP/OPTIMIZER/SEDUMI Optimization
n=size(B,1); m=size(B,2);
Q=sdpvar(n,n,'symmetric');
X=sdpvar(m,m,'symmetric');
Y=sdpvar(m,n,'full');
gama=sdpvar(1);

xk = sdpvar(size(A,1),1);
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
LMIs=[LMIs, [1 xk';xk Q]>=0]; % adds a state restriction
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

    delta = -2*rand + 1; % random number in (-1,1)
    A = [-4 4;-5 0] + delta*[-2 2;-1 4]; % updates the time varying model
    QY = model{x(:,k)}; % updates the Y and Q decision variables
    F{k} = QY{1}/(QY{2}); % is equal to the static gain for k = 1

    u(k) = F{k}*(x(:,k));
    x(:,k+1) = A*x(:,k) + B*u(k);
    y(k+1) = C*x(:,k+1);
    us(k) = Fs*(xs(:,k));
    xs(:,k+1) = A*xs(:,k) + B*us(k);
    ys(k+1) = C*xs(:,k+1);

end

t = 0:ts:tsim;

% figures
figure(1)
subplot(1,2,1)
plot(t,y,'k', t,ys,'k--');
xlabel('Time (sec)')
ylabel('y')
legend('Robust MPC','Static MPC')
subplot(1,2,2)
plot(t(1:end-1),u,'k-',t(1:end-1),us,'k--');
xlabel('Time (sec)')
ylabel('u')






