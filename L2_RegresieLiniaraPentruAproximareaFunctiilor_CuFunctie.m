% regresia liniara pt aproximarea functiilor
close all
load('lab2_09.mat');
n = 15;

plot(id.X, id.Y, 'g')
hold on
plot(val.X, val.Y, 'b')
legend('identificare', 'validare')

% Identificare
figure
plot(id.X, id.Y, 'r') , hold
y_hat_identificare = cauta(n, id.X, id.Y);
plot(id.X, y_hat_identificare, '--black') , title('Identificare');
legend('idY', 'yhatid')

% eroarea datelor de identificare
eroare_identificare = zeros(n, 1);
for i = 1:n
    coeficienti = cauta(i, id.X, id.Y);
    eroare_identificare(i) = eroare(length(id.Y), id.Y, coeficienti);
end

% afisarea erorii pentru fiecare grad 
disp('Erori pentru identificare: ');
disp(eroare_identificare);
% gradul cu cea mai mica eroare pentru valorile de identificare
[min_mse_identificare, grad_optim_identificare] = min(eroare_identificare);
% gradul cu cea mai mica eroare pentru valorile de identificare
disp(['Gradul optim pentru identificare: ', num2str(grad_optim_identificare)]);
% cea mai mica eroare pentru valorile de identificare
disp(['Eroarea minima pentru identificare: ', num2str(min_mse_identificare)]);

% Validare
figure
plot(val.X, val.Y, 'r'), hold on
y_hat_validare = cauta(n, val.X, val.Y);
plot(val.X, y_hat_validare, 'black'), title('Validare');
legend('valY', 'yhatval')

% eroarea datelor de validare
figure
eroare_validare_valori = zeros(n, 1);
for i = 1:n
    coeficienti = cauta(i, val.X, val.Y);
    eroare_validare = eroare(length(val.X), val.Y, coeficienti);
    eroare_validare_valori(i) = eroare_validare;
end
plot(1:n, eroare_validare_valori);
title('MSE pentru validare');
% gradul cu cea mai mica eroare pentru valorile de identificare
[min_mse_validare, grad_optim_validare] = min(eroare_validare_valori); 
% gradul cu cea mai mica eroare pentru valorile de validare
disp(['Gradul optim pentru validare ', num2str(grad_optim_validare)]);
% cea mai mica eroare pentru valorile de validare
disp(['Eroarea minima pentru validare: ', num2str(min_mse_validare)]);

figure
plot(val.X, val.Y, val.X, y_hat_validare);
title('Validare vs Aproximare');

function y_hat = cauta(n, x, y) % creaza y_hat
    phi = [];
    for i = 1:n
        phi = [phi (x(1,1:length(x))') .^ (i-1)];
    end
    theta = phi \ y';
    y_hat = phi * theta;
end

function MSE = eroare(n, y, y_hat) % functie pentru eroarea medie patratica
    e = y - y_hat' ;
    MSE = sum(e(:) .^ 2) / n;
end
