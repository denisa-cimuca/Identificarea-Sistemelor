% regresia liniara pt aproximarea functiilor
close all

load('lab2_09.mat')
%n=15;

plot(id.X, id.Y, 'g')
hold on
plot(val.X, val.Y, 'b')
legend('identificare', 'validare')

MSE__vector_id = [];
MSE_vector_val = [];

for n = 1:22
phi_id = [];
    for i = 1:n
        phi_id = [phi_id (id.X(1, 1:length(id.X))').^(i-1)];
    end
theta = phi_id \ id.Y';
y_hat_id = phi_id*theta;

phi_val = [];
    for i=1:n
        phi_val = [phi_val (val.X(1, 1:length(val.X))') .^ (i-1)];
    end
y_hat_val = phi_val*theta;

N_id = length(id.Y);
suma_id = 0;
    for i=1:N_id
        suma_id = suma_id + (id.Y(i)-y_hat_id(i)).^2;
    end
N_val = length(val.Y);
suma_val = 0;
    for i=1:N_val
        suma_val = suma_val + ((val.Y(i) - y_hat_val(i))).^2;
    end

mse_identificare = suma_id / N_id;
MSE__vector_id(n) = mse_identificare;

mse_validare = suma_val / N_val;
MSE_vector_val(n) = mse_validare;

end

figure
plot(id.X, id.Y)
hold on
plot(id.X, y_hat_id,'r')
title('Aproximare identificare')

figure
plot(val.X, val.Y)
hold on
plot(val.X, y_hat_val, 'r')
title('Aproximare validare')

figure
plot(MSE__vector_id)

figure
plot(MSE_vector_val)

[mse_min_id, index_id] = min(MSE_vector_id)
[mse_min_val, index_val] = min(MSE_vector_val)
