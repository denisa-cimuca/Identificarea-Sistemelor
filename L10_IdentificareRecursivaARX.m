clear all
close all
clc

% u_z = zeros(20,1);
% u_step = 0.3*ones(70,1);
% u_val = [u_z;u_step;u_z];
% u_id = idinput(200,'prbs',[],[-0.8 0.8]);

u_id = idinput(200, 'prbs', [], [-0.8 0.8]);
u_val = [zeros(50,1); 0.3*ones(70,1);zeros(50,1)];

na = 6;
nb = 6;
Ts = 0.01;

motor = DCMRun.start("Ts", 10e-3);

for k = 1:length(u_val)
        motor.wait;
    y_val(k)= motor.step(u_val(k));
end

theta_hat = zeros(na + nb,1);
p_inv = 1000*eye(na + nb,na + nb);
eroare = [];
w = zeros(na+nb);

y_id = [];
for k = 1:length(u_id)
    motor.wait;
    y_id(k) = motor.step(u_id(k));
   
   PHI = zeros(na+nb,1);
    for i = 1:na 
        if k-i <= 0
            PHI(i) = 0;
        else
            PHI(i) = -y_id(k-i);
        end
        
    end

    for i = na + 1:na + nb
        c = i - na;
        if k-c <= 0
            PHI(i) = 0;
        else
            PHI(i) = u_id(k-c);
        end
    end
    if k == 1
        eroare = y_id(k);
    else
    eroare = y_id(k) - PHI'*theta_hat;
    end
    PHI = PHI';
    invP = invP - (invP*(phi')*phi*invP)/(1+phi*invP*(phi'));
    w = p_inv*PHI';  %vector coloana de na+nb
    theta_hat = theta_hat + w*eroare;

    motor.wait();
end
motor.stop();

A = zeros(1, na);
A(1) = 1;
for i = 1 : na 
    A(i+1) = theta_hat(i);
end

B = zeros(1, nb);
B(1) = 0;
for i = 1 : nb
    B(i+1) = theta_hat(i + na);
end

model_complet = idpoly(A,B,[],[],[],0,Ts);
val = iddata(y_val',u_val,Ts);

figure
compare(model_complet,val)