clear all
close all

% N = 200
% u = [zeros(50, 1); idinput(N, 'prbs', [], [-0.7 0.7]); zeros(50,1); 0.4*ones(70,1)] ;
% [vel, alpha, t] = run(u, '3');
% plot(t, vel);

load('date_lab8.mat')

plot(u); title("Intrarea")
figure
plot(vel); title("Iesirea")

u_id = u(1 : 360);
vel_id = vel(1 : 360);

u_VAL = u(450 : 650)';
vel_VAL = vel(450 : 650)';

Ts = 0.05;
data_val = iddata(vel_VAL, u_VAL', Ts);

% parametri
alfa = 0.1; % pragul
Lmax = 27; % limita
delta = 0.001; % gradul de convergenta

nk = 2; % intarzierea
N_u = length(u_id);
% erorile
eroare = zeros(1, N_u);
eroare(1) = 0;

theta = [-1,1.99];
theta_urmator = [1,1];
% b = theta(2);
% f = theta(1);

% derivatele
def = zeros(N_u, 1);
deb = zeros(N_u, 1);
deriv_e = zeros(N_u, 2);
% hessianul
Hessian = zeros(2, 2);
% gradientul
gradient = zeros(2, 1);

l = 0;
while l<Lmax % || norm(thata - theta_urmator)>delta
    for k = 1:nk 
        eroare(k) = vel_id(k);
        deriv_e(1, k)=0;
        deriv_e(2, k) = 0;
    end
    for k = nk+1:N_u
        eroare(k) = vel_id(k) + vel_id(k - 1)*theta(1) - eroare(k - 1)*theta(1) - u_id(k - nk)*theta(2);
        def(k)=-def(k - 1)*theta(1) + vel_id(k - 1) - eroare(k - 1);
        deb(k) = -deb(k-1)*theta(1) - u_id(k-nk);
    end
deriv_e = [def deb]';

    var = zeros(2, 1);
    for k=1:N_u
        gradient = gradient + (eroare(k)*[def(k) deb(k)]');
        var = [deriv_e(1,k);deriv_e(2,k)];
        Hessian = Hessian + var*var';
    end
    gradient = (2/(N_u - nk))*gradient;
    Hessian = (2/(N_u - nk))*Hessian;

    theta_urmator= theta - alfa*inv(Hessian)*gradient;
    l = l + 1;
    norma = norm(theta - theta_urmator);
% if((morma>=delta)||l==lmax)
%     break;
% end
theta = theta_urmator;
end

f = theta_urmator(1);
b = theta_urmator(2);

figure
model = idpoly(1,[0,b],1,1,[1,f],0,Ts);
compare(model, data_val)