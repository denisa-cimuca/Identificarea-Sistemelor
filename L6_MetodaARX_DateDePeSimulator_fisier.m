clear all
load('lab6_8.mat')

subplot(211)
plot(id)
subplot(212)
plot(val)
title('Intrarea si Iesirea')

u_id = id.InputData;
u_val = val.InputData;
vel_id = id.OutputData;
vel_val = val.OutputData;

n = length(u_id);
na = 10;
nb = 10;
phi = [];
for k = 1:n
    for i = 1:na
        if k>i
            phi(k,i) = -vel_id(k-i);
        else
            phi(k,i) = 0;
        end
    end
    for j = 1:nb
        if k>j
            phi(k,j+na) = u_id(k-j);
        else
            phi(k,j+na) = 0;
        end
    end
end

theta = phi \ vel_id;
%yhat_id = phi*theta;

% figure
% plot(vel_id)
% hold on
% plot(yhat_id)
% title('Predictie identificare')


n_val = length(u_val);
% phi_val = [];
% for k = 1:n_val
%     for i = 1:na
%         if k>i
%             phi_val(k,i) = -vel_val(k-i);
%         else
%             phi_val(k,i) = 0;
%         end
%     end
%     for j = 1:nb
%         if k>j
%             phi_val(k,j+na) = u_val(k-j);
%         else
%             phi_val(k,j+na) = 0;
%         end
%     end
% end
% 
% yhat_val = phi_val*theta;
% 
% figure
% plot(vel_val)
% hold on
% plot(yhat_val)
% title('Predictie validare')


% % Simularea
% yhat_id_simulat = zeros(length(u_id), 1);
% for k =1:length(u_id)
%     for i = 1:na
%             if(k>i) 
%                 yhat_id_simulat(k) = yhat_id_simulat(k) + (-yhat_id_simulat(k-i)*theta(i)); 
%             end
%         end
%         for j = 1:nb
%             if(k>j) 
%                 yhat_id_simulat(k) = yhat_id_simulat(k) + u_id(k-j)*theta(na+j);
%             end
%         end
%         %yhat_id_simulat(k) = suma;
% end
% 
% figure
% plot(vel_id)
% hold on
% plot(yhat_id_simulat)
% title('Simulare identificare')


yhat_val_simulat = zeros(length(u_val),1);
    for k = 1:length(u_val)
        for i = 1:na
            if(k>i) 
                yhat_val_simulat(k) = yhat_val_simulat(k) + (-yhat_val_simulat(k-i)*theta(i)); 
            end
        end
        for j = 1:nb
            if(k>j) 
                yhat_val_simulat(k) = yhat_val_simulat(k) + u_val(k-j)*theta(na+j);
            end
        end
        %yhat_val_simulat(k) = suma;

    end

figure
plot(vel_val)
hold on
plot(yhat_val_simulat)
title('Simulare validare')


%e = vel_val-yhat_val_simulat;
MSE = sum((yhat_val_simulat-vel_val).^2) * 1/length(yhat_val_simulat)
%min_mse=min(MSE)


figure
na = 1:20;
nb = 1:20;
nk = 1;

NN = struc(na,nb,nk);

V = arxstruc(id,val,NN);
N = selstruc(V, 0);

model = arx(id,N);

compare(val,model);