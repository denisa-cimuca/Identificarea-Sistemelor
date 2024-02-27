% analiza de corelatia
load('lab5_8.mat')

u_id = id.InputData;
y_id = id.OutputData;

u_val = val.InputData;
y_val = val.OutputData;

plot(detrend(id))

MSE_vector_id =[];
MSE_vector_val = [];

N = length(u_id);

ru = [];
for tau = 0:N-1
    suma = 0;
    for k = 1:N-tau
        suma = suma + u_id(k+tau)*u_id(k);
    end
    ru(tau+1) = suma/N;
end

ryu = [];
for tau = 0:N-1
    suma = 0;
    for k =1:N-tau
        suma = suma + y_id(k+tau)*u_id(k);
    end
    ryu(tau+1) = suma/N;
end

M = 49;
RU = [];
for i = 1:N
    for j = 1:M
        RU(i,j) = ru(abs(i-j)+1);
    end
end

h = RU\ryu';

y_hat_id = conv(h, u_id);
y_hat_id = y_hat_id(1:length(y_id));

y_hat_val = conv(h, u_val);
y_hat_val = y_hat_val(1:length(y_val));

figure
plot(tid, y_id)
hold on
plot(tid, y_hat_id,'g')

figure
plot(tval, y_val)
hold on
plot(tval, y_hat_val, 'g')

e = y_val-y_hat_val;
MSE = sum(e.^2)/length(e)


figure
plot(imp)
hold on
plot(h,'g')