% ARX 
% model 1 cu validare 1 : pentru m=3
load('date_lab7.mat') 

plot(u)
figure
plot(vel)

u_id_m3 = u(1:257);
vel_id_m3 = vel(1:257);

u_val = u(547:end);
vel_val = vel(547:end);

na_m3 = 3;
nb_m3 = 5;

% Predictie
Nid_m3 = length(u_id_m3); % lungimea intrarii identificate

ID_phi_m3 = [];
for tau_id = 1 : Nid_m3
    for i = 1 : na_m3
        if tau_id > i
            ID_phi_m3(tau_id, i) = -vel_id_m3(tau_id - i);
        else
            ID_phi_m3(tau_id , i) = 0;
        end
    end
    for j = 1 : nb_m3
        if tau_id > j
            ID_phi_m3(tau_id, j+na_m3) = u_id_m3(tau_id - j);
        else
            ID_phi_m3(tau_id, j+na_m3) = 0;
        end
    end
end

theta_m3 = ID_phi_m3 \ vel_id_m3';
yhat_id_pred_m3 = ID_phi_m3 * theta_m3;

% figure
% hold on
% plot(vel_id, 'b')
% plot(yhat_id_pred_m3, 'r')
% title('predictie identificare')

Nval_m3 = length(u_val); % lungimea intrarii de validare

VAL_phi_m3 = [];
for tau_val = 1 : Nval_m3
    for i = 1 : na_m3
        if tau_val > i
            VAL_phi_m3(tau_val, i) = -vel_val(tau_val - i);
        else
            VAL_phi_m3(tau_val, i) = 0;
        end
    end
    for j = 1 : nb_m3
        if tau_val > j
            VAL_phi_m3(tau_val, j+na_m3) = u_val(tau_val - j);
        else
            VAL_phi_m3(tau_val, j+na_m3) = 0;
        end
    end
end

yhat_val_pred_m3 = VAL_phi_m3 * theta_m3; % iesirea prin predicite


%Simulare
yhat_val_Simulat_m3 = zeros(length(u_val), 1); % iesirea prin simulare
for tau = 1 : Nval_m3
    for i = 1 : na_m3
        if tau > i
            yhat_val_Simulat_m3(tau) = yhat_val_Simulat_m3(tau) - yhat_val_Simulat_m3(tau - i)*theta_m3(i);
        end
    end
    for j = 1 : nb_m3
        if tau > j 
            yhat_val_Simulat_m3(tau) = yhat_val_Simulat_m3(tau) + u_val(tau - j)*theta_m3(j + na_m3);
        end
    end
end

% Predictia
figure; title('Predictie validare')
hold on
plot(vel_val, 'b')
plot(yhat_val_pred_m3, 'r')

% Simularea
figure; title('Simulare validare')
hold on
plot(vel_val, 'b')
plot(yhat_val_Simulat_m3, 'r')

MSE_m3 = sum((yhat_val_Simulat_m3 - vel_val') .^ 2) * 1 / length(yhat_val_Simulat_m3)


%% model 2 cu validare 2 : pentru pentru m=10
% ARX 
load('date_lab7.mat')
%load('data1_lab7.mat')

plot(u)
figure
plot(vel)

u_id_m10 = u(300:500);
vel_id_m10 = vel(300:500);

u_val = u(547:end);
vel_val = vel(547:end);

na_m10 = 10;
nb_m10 = 10;

% Predictie
Nid_m10 = length(u_id_m10); % lungimea intrarii identificate

ID_phi_m10 = [];
for tau_id = 1 : Nid_m10
    for i = 1 : na_m10
        if tau_id > i
            ID_phi_m10(tau_id, i) = -vel_id_m10(tau_id - i);
        else
            ID_phi_m10(tau_id , i) = 0;
        end
    end
    for j = 1 : nb_m10
        if tau_id > j
            ID_phi_m10(tau_id, j+na_m10) = u_id_m10(tau_id - j);
        else
            ID_phi_m10(tau_id, j+na_m10) = 0;
        end
    end
end

theta_m10 = ID_phi_m10 \ vel_id_m10';
yhat_id_pred_m10 = ID_phi_m10 * theta_m10;

% figure
% hold on
% plot(vel_id, 'b')
% plot(yhat_id_pred_m10, 'r')
% title('predictie identificare')

Nval_m10 = length(u_val); % lungimea intrarii de validare

VAL_phi_m10 = [];
for tau_val = 1 : Nval_m10
    for i = 1 : na_m10
        if tau_val > i
            VAL_phi_m10(tau_val, i) = -vel_val(tau_val - i);
        else
            VAL_phi_m10(tau_val, i) = 0;
        end
    end
    for j = 1 : nb_m10
        if tau_val > j
            VAL_phi_m10(tau_val, j+na_m10) = u_val(tau_val - j);
        else
            VAL_phi_m10(tau_val, j+na_m10) = 0;
        end
    end
end

yhat_val_pred_m10 = VAL_phi_m10 * theta_m10; % iesirea prin predicite


%Simulare
yhat_val_Simulat_m10 = zeros(length(u_val), 1); % iesirea prin simulare
for tau = 1 : Nval_m10
    for i = 1 : na_m10
        if tau > i
            yhat_val_Simulat_m10(tau) = yhat_val_Simulat_m10(tau) - yhat_val_Simulat_m10(tau - i)*theta_m10(i);
        end
    end
    for j = 1 : nb_m10
        if tau > j 
            yhat_val_Simulat_m10(tau) = yhat_val_Simulat_m10(tau) + u_val(tau - j)*theta_m10(j + na_m10);
        end
    end
end

% Predictia
figure; title('Predictie validare')
hold on
plot(vel_val, 'b')
plot(yhat_val_pred_m10, 'r')

% Simularea
figure; title('Simulare validare')
hold on
plot(vel_val, 'b')
plot(yhat_val_Simulat_m10, 'r')

MSE_m10 = sum((yhat_val_Simulat_m10 - vel_val') .^ 2) * 1 / length(yhat_val_Simulat_m10)

%% u total
u = [zeros(1,50), u3, zeros(1,50), u10, zeros(1,50), 0.4*ones(1,70)];

% [vel, alpha, t] = run(u, '4');
% plot(t,vel);

plot(u)
%% determinarea intrarilor 
% dimensiunea matricei A3, m=3
m3=3;
% generarea matricei A3
A3=zeros(m3,m3);
linie_m3 = [];

if  m3==3
 linie_m3= [1 0 1];
elseif m3==4
    linie_m3=[1 0 0 1];
elseif m3==5
    linie_m3=[0 1 0 0 1];
elseif m3==6
    linie_m3=[1 0 0 0 0 1];
elseif m3==7
    linie_m3=[1 0 0 0 0 0 1];
elseif m3==8
    linie_m3=[1 1 0 0 0 0 0 1];
elseif m3==9
    linie_m3=[0 0 0 1 0 0 0 0 1];
elseif m3==10
    linie_m3=[0 0 1 0 0 0 0 0 0 1];
end

A3(1,:)=linie_m3;

for i = 2:m3
    A3(i,i-1)=1;
end

%  dimensiunea matricei A10, m=10
m10=10;
% generarea matricei A10
A10=zeros(m10,m10);
linie_m10=[];

if  m10==3
 linie_m10= [1 0 1];
elseif m10==4
    linie_m10=[1 0 0 1];
elseif m10==5
    linie_m10=[0 1 0 0 1];
elseif m10==6
    linie_m10=[1 0 0 0 0 1];
elseif m10==7
    linie_m10=[1 0 0 0 0 0 1];
elseif m10==8
    linie_m10=[1 1 0 0 0 0 0 1];
elseif m10==9
    linie_m10=[0 0 0 1 0 0 0 0 1];
elseif m10==10
    linie_m10=[0 0 1 0 0 0 0 0 0 1];
end

A10(1, :)=linie_m10;

for i = 2 : m10
    A10(i, i-1)=1;
end

% numar de esantioane
N = 200;
% limite
a = -0.7;
b = 0.7;

% SPAB pentru m=3
u3  = creareSPAB(N, m3, a, b, A3); %% PE = persistenta a excitatiei
figure, plot(u3);
% SPAB pentru m=10
u10  = creareSPAB(N, m10, a, b, A10); %% PE = persistenta a excitatiei
figure, plot(u10);

function u = creareSPAB(N, m, a, b, A)
    %P = 2^m-1; perioada = P;
    %X(1,:)=ones(1,m);
    Q = zeros(N, m);
    Q(1, :) = 1;
    for k = 1:N-1        
            Q(k+1, :) = (mod(A * Q(k, :)', 2));
        u(k) = Q(k, m);
    end
    % gererarea intrarii u
    u(m) = Q(N, m);
    for i = 1 : N-1
        u(i) = a + (b - a)*u(i); 
    end
end

%% CU FUNCTIA ARX 

% arx pe m3 cu functia arx
load('date_lab7.mat')

subplot(211),plot(u)
subplot(212),plot(vel)

u_id_m3 = u(50:272);
vel_id_m3 = vel(50:272);
Ts = 0.05;
data_id = iddata(vel_id_m3', u_id_m3',Ts);

na = 3;
nb = 5;
nk=1;
model_identificare = arx(data_id, [na,nb,nk]);

u_val = u(540:end);
vel_val = vel(540:end);
data_val = iddata(vel_val', u_val',Ts);

figure
compare(model_identificare,data_val)

% arx pe m10 cu functia arx
load('date_lab7.mat')

subplot(211),plot(u)
subplot(212),plot(vel)

u_id_m10 = u(299:500);
vel_id_m10 = vel(299:500);
Ts = 0.05;
data_id_m10 = iddata(vel_id_m10', u_id_m10',Ts);

na = 10;
nb = 10;
nk=1;
model_identificare_m10 = arx(data_id_m10, [na,nb,nk]);

u_val = u(540:end);
vel_val = vel(540:end);
data_val = iddata(vel_val', u_val',Ts);

figure
compare(model_identificare_m10,data_val)