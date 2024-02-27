close all
clear all

%extragere date de pe motor
%N = 200;
%u = [zeros(50, 1); idinput(N, 'prbs', [], [-0.7 0.7]); zeros(50,1); 0.4*ones(70,1)] ;
% [vel, alpha, t] = run(u, '3');
% plot(t, vel);

load('date_lab6') % de pe motor

subplot(211)
plot(u)
subplot(212)
plot(vel)
title('Intrarea si Iesirea')

u_id = u(51 : 250);
vel_id = vel(51 : 250);

u_val = u(351 : 550);
vel_val = vel(351 : 550);

na = 10;
nb = 10;

% Predictie
Nid = length(u_id); % lungimea intrarii identificate

ID_phi = [];
for tau_id = 1 : Nid
    for i = 1 : na
        if tau_id > i
            ID_phi(tau_id, i) = -vel_id(tau_id - i);
        else
            ID_phi(tau_id , i) = 0;
        end
    end
    for j = 1 : nb
        if tau_id > j
            ID_phi(tau_id, j+na) = u_id(tau_id - j);
        else
            ID_phi(tau_id, j+na) = 0;
        end
    end
end

theta = ID_phi \ vel_id';
yhat_id_pred = ID_phi * theta;

% figure
% hold on
% plot(vel_id, 'b')
% plot(yhat_id_pred, 'r')
% title('predictie identificare')

Nval = length(u_val); % lungimea intrarii de validare

VAL_phi = [];
for tau_val = 1 : Nval
    for i = 1 : na
        if tau_val > i
            VAL_phi(tau_val, i) = -vel_val(tau_val - i);
        else
            VAL_phi(tau_val, i) = 0;
        end
    end
    for j = 1 : nb
        if tau_val > j
            VAL_phi(tau_val, j+na) = u_val(tau_val - j);
        else
            VAL_phi(tau_val, j+na) = 0;
        end
    end
end

yhat_val_pred = VAL_phi * theta; % iesirea prin predicite


%Simulare
yhat_val_Simulat = zeros(length(u_val), 1); % iesirea prin simulare
for tau = 1 : Nval
    for i = 1 : na
        if tau > i
            yhat_val_Simulat(tau) = yhat_val_Simulat(tau) - yhat_val_Simulat(tau - i)*theta(i);
        end
    end
    for j = 1 : nb
        if tau > j 
            yhat_val_Simulat(tau) = yhat_val_Simulat(tau) + u_val(tau - j)*theta(j + na);
        end
    end
end

% Predictia
figure; title('Predictie validare')
hold on
plot(vel_val, 'b')
plot(yhat_val_pred, 'r')

% Simularea
figure; title('Simulare validare')
hold on
plot(vel_val, 'b')
plot(yhat_val_Simulat, 'r')

MSE = sum((yhat_val_Simulat - vel_val') .^ 2) * 1 / length(yhat_val_Simulat)