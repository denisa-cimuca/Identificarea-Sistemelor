% Ordinul I
clear all
% analiza raspunsurilor la impuls
load('lab4_order1_8.mat')

u = data.u;
y = data.y;

plot(t,u, 'g')
hold on
plot(t,y, 'r')

yss = sum(y(120:130)) / length(y(120:130))
y0 = yss
ymax = max(y)
uss = u(1)
u0 = uss
umax = max(u)
k = yss/uss;
y_descreste = y0 + 0.368*(ymax-y0)
%[ya, t2] = min(abs(y(1:100)-y_descreste))
t1 = t(32);  % timpul unde este ymax = de unde pleaca t in sus
t2 = t(45);  % timpul unde este y_descreste
T = t2 - t1;

A = -1/T;
B = k/T;
C = 1;
D = 0;
Htf = ss(A,B,C,D);

y_id = lsim(Htf, u(1:130), t(1:130), y0);
y_val = lsim(Htf, u(131:330), t(131:330), y0);

figure
plot(t, y, 'b')
hold on
plot(t(1:130), y_id, 'g')
hold on
plot(t(131:330), y_val, 'r')

e = y(131:330) - y_val;
MSE = sum(e.^2)/length(e)

%% Ordinul II
% analiza raspunsurilor la impuls
load('lab4_order2_7.mat')
 
u = data.u;
y = data.y;

plot(t,u)
hold on
plot(t,y)
legend('Intrarea', 'Iesirea')

Ts = t(9)-t(8);
yss = sum(y(120:130),1)/length(y(120:130))
y0 = yss
ymax = max(y)
uss = u(1)
u0 = uss
umax = max(u)
k = yss/uss

t_plecare = 1.5; % adica t(31)
t_intersectie_2 = 2.4; % adica t(49) cu axa de pornire a lui y
t_intersectie_3 = 3.25; % adica t(66)

% varfurile iesirii
t1 = 1.85; 
t2 = 2.65;
t3 = 3.75;

T0 = t3 - t1;

% k-urile pentru formula ariilor Ap, Am: Ap = Ts * sum(y(k) - y0)
indice1 = 31;
indice2 = 49;
indice3 = 66;

Aplus = Ts * sum(y(indice1:indice2)-y0);
Aminus = Ts * sum(y0 - y(indice2:indice3));
M = Aminus/Aplus;

zeta = (log(1 / M)) / sqrt(pi ^ 2 + (log(M)) ^ 2);
wn = 2*pi/(T0*sqrt(1-zeta^2));

A = [0 1; -wn^2 -2*zeta*wn];
B = [0; k*wn^2];
C = [1 0];
D = 0;

Htf = ss(A,B,C,D);

y_id = lsim(Htf, u(1:130), t(1:130), [y0 0]);
y_val = lsim(Htf, u(131:330), t(131:330), [y0 0]);

figure
plot(t, y, 'b')
hold on
plot(t(1:130), y_id, 'g')
hold on
plot(t(131:330), y_val, 'r')

e = y(131:330) - y_val;
MSE = sum(e.^2) / length(e)


% aria sub raspunsul sistemului
aria_raspuns = trapz(t, y);
% aria de suprareglaj
deriv = y-yss;
aria_suprareglaj = trapz(t, deriv)
