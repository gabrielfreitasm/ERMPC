clear, clc, close all
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

%% Model information
% y(t) = (G1d)u1(t) + (G2d)u2(t)
Ts = 0.02;
b1b = [0 0.0318 0.1187];
a1 = [1 -1.8195 0.8239];
b2b = [0 -0.0239 -0.0673];
a2 = [1 -1.7215 0.7325];
G1d = tf(b1b,a1,Ts);
G2d = tf(b2b,a2,Ts);
Gd = [G1d G2d];
[A,B,C,D] = ssdata(Gd);

%% Polytope vertices

A1 = A - 0.2*A;
A2 = A + 0.2*A;

%% Control

% Project restrictions
umax = 500;

% Parameters Qc Rc
Qc = 1; Rc = 2e-3;

% ONLINE Model
% YALMIP/OPTIMIZER/SEDUMI Optimization
n=size(B,1); m=size(B,2);
Q=sdpvar(n,n,'symmetric');
X=sdpvar(m,m,'symmetric');
Y=sdpvar(m,n,'full');
gama=sdpvar(1);

xk = sdpvar(4,1);
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
x = [1 0 0 0]';
xs = [1 0 0 0]';
QY = model{x}; % returns Y and Q for the initial condition x
Fs = QY{1}/(QY{2}); % static model

% Simulation parameters
tsim = 10;
num = ceil(tsim/Ts);
delta = zeros(1,num);
u = zeros(m,num);
us = zeros(m,num);
F = cell(1,num); % dynamic gain - calculated in each iteraction

for k = 1:num

    delta(k) = -0.1*rand + 0.05; % random number in (-0.2,0.2)
    A = A + delta(k)*A; % updates the time varying model
    QY = model{x(:,k)}; % updates the Y and Q decision variables
    F{k} = QY{1}/(QY{2}); % is equal to the static gain for k = 1

    u(:,1) = F{k}*(x(:,k));
    x(:,k+1) = A*x(:,k) + B*u(:,1);
    y(k+1) = C*x(:,k+1);
    us(:,1) = Fs*(xs(:,k));
    xs(:,k+1) = A*xs(:,k) + B*us(:,1);
    ys(k+1) = C*x(:,k+1);

end

t = 0:Ts:tsim;

% figures
figure(1)
subplot(1,2,1)
plot(t,y,'k', t,ys,'k--'); % output
xlabel('Time (sec)','FontSize',15)
ylabel('Output','FontSize',15)
legend('Robust MPC','Static MPC','FontSize',12)
subplot(1,2,2)
plot(t(1:end-1),u(1,:),'b-',t(1:end-1),us(1,:),'b--');
hold on
plot(t(1:end-1),u(1,:),'k-',t(1:end-1),us(1,:),'k--');
xlabel('Time (sec)','FontSize',15)
ylabel('u (volts)','FontSize',15)






