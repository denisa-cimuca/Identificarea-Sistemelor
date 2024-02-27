clear all
close all
load('date_lab9.mat')

% N=300;
% u = [zeros(50, 1); idinput(N, 'prbs', [], [-0.8 0.8]); zeros(100,1); 0.3*ones(70,1)] ;
% [vel, alpha, t] = run(u, '3');
% plot(t, vel);

plot(u); title("Intrarea")
figure
plot(vel); title("Iesirea")

Ts = 0.05;
na = 9;
nb = 9;
nk = 1;

u_id = u(1 : 360);
vel_id = vel(1 : 360)';
data_id = iddata(vel_id, u_id, Ts);
model_identificare = arx(data_id, [na, nb, nk]);

u_val = u(450 : 650);
vel_val = vel(450 : 650)';
data_val = iddata(vel_val, u_val, Ts);
model_validare = arx(data_val, [na, nb, nk]);

% Simulare
yhat_simulat = lsim(model_identificare, u_id);

Phi_simulare = [];
for tau = 1 : length(yhat_simulat)
    for i = 1 : na
        if tau > i
            Phi_simulare(tau, i) = -yhat_simulat(tau - i);
        else
            Phi_simulare(tau, i) = 0;
        end
    end
    for j = 1 : nb
        if tau > j
            Phi_simulare(tau, j + na) = u_id(tau - j);
        else
            Phi_simulare(tau, j + na) = 0;
        end
    end
end

theta = Phi_simulare \ yhat_simulat;

A = zeros(1, na+1);
A(1) = 1;
for i = 1 : na 
    A(i+1) = theta(i);
end

B = zeros(1, nb+1);
B(1) = 0;
for i = 1 : nb
    B(i+1) = theta(i + na);
end

ivmodel = idpoly(A, B, 1, 1, 1, 0 ,Ts);
yhat = lsim(ivmodel, u_val);
% compare(ivmodel, data_val)

figure
hold on
plot(yhat, 'r');
plot(vel_val, 'b')