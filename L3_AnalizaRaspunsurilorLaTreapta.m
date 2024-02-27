%Ordin I
clear all
load('lab3_order1_8.mat')

u = data.u;
y = data.y;

subplot(211)
plot(t,u), title('Intrarea')
subplot(212)
plot(t,y), title('iesirea')

% k = y stationar ; 0.632*k gaseste valoarea timpului la care iesirea atinge 0.632 din valoarea stationara

yss = sum(y(90:100))/length(90:100)
%media(90:100)
% length(90:100)=11 este timpul la care k este de 0.632 => x=11, y=8.99

y0 = y(1)  % punctul de plecare al iesirii
uss = u(1) % punctul de plecare al intrarii
u0 = 0 
k = (yss-y0)/(uss-u0);
yT = y0 + 0.632 * (yss-y0); % noua iesire; valoarea iesirii unde y atinge 0.632 din yss

[yT, indice] = min(abs(y(1:100) - yT));
T = t(indice); % valoarea lui T unde iesirea atinge 0.632 din valoarea stationara

% Modelul sistemului 
H = tf(k, [T 1])
y_id = lsim(H, u(1:100), t(1:100));
y_val = lsim(H, u(201:500), t(201:500));

e = y(201:500) - y_val;
MSE = sum(e.^2)/ length(e)

figure
plot(t, y, 'b')
hold on
plot(t(1:100), y_id, 'g')
hold on
plot(t(201:500), y_val,'r')

%% Ordinul II
% analiza raspunsului la treapta
load('lab3_order2_8.mat')

u = data.u;
y = data.y;

subplot(211)
plot(t,u)
subplot(212)
plot(t,y)

yss = sum(y(90:100)) / length(90:100)
y0 = y(1)
uss = u(1)
u0 = 0
K = (yss - y0) / (uss - u0);
yT = y0 + 0.632 * (yss - y0) % valoarea iesirii unde atinge 0.632 din yss
[ya, indice1] = min(abs(y(1:100) - yT))
T = t(indice1)

[yFirstMax, t1] = max(y(1:100))

suprareglaj = (yFirstMax - yss) /yss;
zeta = (log(1 / suprareglaj)) / sqrt(pi ^ 2 + log(suprareglaj) ^ 2);

yFirstMin = yss * (1 - suprareglaj ^ 2)
[yb, t2] = min(abs(y(1:100) - yFirstMin))

T0 = 2 * (t(t2) - t(t1))
wn = (2 * pi) / (T0 * sqrt(1 - zeta ^ 2));

H = tf(K*wn^2, [1 2*zeta*wn wn^2])

y_id = lsim(H, u(1:100), t(1:100));
y_val = lsim(H, u(201:500), t(201:500));

figure
plot(t, y, 'b');
hold on
plot(t(1:100), y_id, 'g');
hold on
plot(t(201:500), y_val, 'r');
title('Sistem de ordin II');

e = y(201:500) - y_val;
MSE2 = sum(e .^ 2) / length(e); 