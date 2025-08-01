clc
clear
warning off
tic
%Initialize Variables
UB=1e5;
LB=1e-5;
k=1;



%1st Iteration
%Initialize forecasts PV and variable load power 
price =[0.47,0.47,0.47,0.47,0.47,0.47,0.47,0.9,1.35,1.35,1.35,0.9,0.9,0.9,0.9,0.9,0.9,0.9,1.35,1.35,1.35,1.35,1.35,0.47];

PPV_f=[0,0,0,0,0,0,0,146,383,711,1070,1120,1246,976,864,928,712,477,279,0,0,0,0,0];


PL_f=[365,327,326,284,275,420,535,546,573,716,763,954,774,699,696,667,768,745,815,882,903,676,505,419];


%Estimated variable load 
PDR_exp=[80,70,60,50,70,70,90,100,120,150,170,200,140,100,100,120,140,150,190,200,200,190,100,80];


%% Parameters
n_days = 100; % 100 days of data
n_hours = 24;

% 1. Electricity Price (¥/kWh)
price_base = [0.47,0.47,0.47,0.47,0.47,0.47,0.47,0.9,1.35,1.35,1.35,0.9,0.9,0.9,0.9,0.9,0.9,0.9,1.35,1.35,1.35,1.35,1.35,0.47];
price = zeros(1, n_days * n_hours);

for day = 1:n_days
    offset = (day-1)*n_hours + 1 : day*n_hours;
    daily_price = price_base;

    % Weekend discount (10% reduction on Saturdays/Sundays)
    if mod(day, 7) == 6 || mod(day, 7) == 0
        daily_price(8:23) = daily_price(8:23) * 0.9; 
    end

    % Random fluctuation (±5% for mid/peak hours)
    noise = 1 + 0.05*(rand(1, 24) - 0.5); 
    daily_price = daily_price .* noise;

    price(offset) = daily_price;
end

% 2. Variable Load (PDR_exp, kW)
PDR_exp_base = [80,70,60,50,70,70,90,100,120,150,170,200,140,100,100,120,140,150,190,200,200,190,100,80];
PDR_exp = zeros(1, n_days * n_hours);

for day = 1:n_days
    offset = (day-1)*n_hours + 1 : day*n_hours;
    daily_PDR = PDR_exp_base;

    % Weekend reduction (15% lower demand)
    if mod(day, 7) == 6 || mod(day, 7) == 0
        daily_PDR = daily_PDR * 0.85;
    end

    % Daily scaling (±10%) and hourly noise (±5%)
    daily_scale = 1 + 0.1*(rand() - 0.5);
    hourly_noise = 1 + 0.05*(rand(1, 24) - 0.5);
    daily_PDR = daily_PDR .* hourly_noise * daily_scale;

    PDR_exp(offset) = daily_PDR;
end

% 3. Load PV,WT,Variable Load Data from Excel

%Parameters
n_5min = 288;   % Fixed data points per day (every 5 minutes)
n_days = 100;     % Number of days to process
sampled_indeces = 1:12:n_5min;


% Load datasets
PL_f_data = xlsread('Load Power Data Over 365 Days.xlsx', 'Sheet1');
PPV_f_data = xlsread('PV Power Data Over 365 Days.xlsx', 'Sheet1'); 
% wt_data = xlsread('WT Power Data Over 365 Days.xlsx', 'Sheet1');


PL_f = PL_f_data(1:n_days,sampled_indeces)';  
PPV_f = PPV_f_data(1:n_days,sampled_indeces)';  
% wt_data = wt_data(1:n_days,sampled_indeces); 


% PL_f = reshape(PL_f',n_hours*n_days, 1);  
% PPV_f= reshape(PPV_f', n_hours*n_days, 1);    
% PL_f_100D = reshape(wt_data,n_hours*n_days, 1); 


% Reshape into 2400x1 vectors
price_100D = price(:);
PDR_exp_100D = PDR_exp(:);
PL_f_100D = 0.25*PL_f(:);
PPV_f_100D = 0.5*PPV_f(:);

disp(size(price_100D));
disp(size(PDR_exp_100D));
disp(size(PL_f_100D));
disp(size(PPV_f_100D));

Cday = [];

%% Optimization Loop

for d = 1:100
    fprintf('%d th day Optimization',d);
    PDR_exp = (PDR_exp_100D(24*(d-1)+1:24*d))';
    price = (price_100D(24*(d-1)+1:24*d))';
    PPV_f = (PPV_f_100D(24*(d-1)+1:24*d))';
    PL_f = (PL_f_100D(24*(d-1)+1:24*d))';

% disp(PPV_f);
% disp(PL_f);
uncert=[PPV_f,PL_f]';




[LB,x,y]=Primal_Problem(k,uncert,price,PDR_exp);
[UB,uncert_SP] = Dual_Problem(k,x,price,PDR_exp,PPV_f,PL_f);
% [UB,u_SP] = SP(k,x,price,Ppv_f,PL_f);
p(1)= UB - LB;
fprintf('%d Iteration\n',k);  
fprintf('Upper Bound：%f\n',UB);
fprintf('Lower Bound：%f\n',LB);
fprintf('Difference ：%f\n',p(k));


% 2nd and more iterations
for k=1:8

    % uncert= uncert_SP;
    uncert=[uncert,uncert_SP];
    [LB,x,y_MP]=Primal_Problem(k+1,uncert_SP,price,PDR_exp); 
    [UB_SP,uncert_SP] = Dual_Problem(k+1,x,price,PDR_exp,PPV_f,PL_f);
    UB=min(UB,UB_SP);
    p(k+1) = UB-LB;
    
fprintf('%d Iteration\n',k+1);  
fprintf('Upper Bound：%f\n',UB);
fprintf('Lower Bound：%f\n',LB);
fprintf('Difference ：%f\n',p(k+1));
    if(p(k+1) < 5)
        break;
    end

k=k+1;
end

Cday(d) = LB;

%% 数据处理
result_y=double(y_MP);
result_u=double(uncert_SP);
result_Ppv(24*(d-1)+1:24*d)=result_u(1:24,1);
result_PL(24*(d-1)+1:24*d)=result_u(25:48,1);

result_PG(24*(d-1)+1:24*d)=result_y(1:24,1);
result_Pch(24*(d-1)+1:24*d)=result_y(25:48,1);
result_Pdis(24*(d-1)+1:24*d)=result_y(49:72,1);
result_PDR(24*(d-1)+1:24*d)=result_y(73:96,1);
result_Pbuy(24*(d-1)+1:24*d)=result_y(145:168,1);
result_Psell(24*(d-1)+1:24*d)=result_y(169:192,1);
end

disp(sum(result_Pch(24*99+1:24*100)));
disp(sum(result_Pdis(24*99+1:24*100)));
disp(sum(result_Pbuy(24*99+1:24*100)));
disp(sum(result_Psell(24*99+1:24*100)));
disp(sum(result_Ppv(24*99+1:24*100)));
disp(sum(result_PL(24*99+1:24*100)));
disp(sum(result_Ppv));
disp(sum(result_PL));
disp(sum(result_PG));

%% 画图
figure(1)
subplot(2,1,1)

bar(-result_Pch,0.75,'b')
hold on
bar(result_Pdis,0.75,'g')
xlim([1 2400])
ylim([-600 600])
grid
legend('充电功率','放电功率');
xlabel('时间/h')
ylabel('功率/kw')
subplot(2,1,2)
bar(-result_Pbuy,0.75,'b')
grid
hold on
bar(result_Psell,0.75,'g')
xlim([1 2400])
ylim([-1500 1500])
legend('购电功率','售电功率');
xlabel('时间/h')
ylabel('功率/kw')

 figure(4);
 grid
 plot(p(1:k))
 xlim([0 4])
 xlabel('迭代次数')
 ylabel('UB-LB')
 title('运行曲线')

figure(2)
grid
plot(PPV_f_100D','k','linewidth',1)
hold on
plot(result_Ppv,'r.--','linewidth',1)
hold on
plot(1.15*result_Ppv,'g.--','linewidth',1)
hold on
plot(0.85*result_Ppv,'g.--','linewidth',1)
xlim([1 2400])
ylim([0 1600])
legend('光伏实际出力','光伏预测出力','光伏区间出力上限','光伏区间出力下限');
xlabel('时间/h')
ylabel('功率/kw')

figure(3)
grid
plot(PL_f_100D','k','linewidth',1)
hold on
plot(result_PL,'r.--','linewidth',1)
hold on
plot(1.1*result_PL,'g.--','linewidth',1)
hold on
plot(0.9*result_PL,'g.--','linewidth',1)
xlim([1 2400])
ylim([0 1200])
legend('负荷实际出力','负荷预测出力','负荷区间出力上限','负荷区间出力下限');
xlabel('时间/h')
ylabel('功率/kw')

figure(5)

bar([result_PDR,PDR_exp_100D']);
grid
xlim([1 2400])
ylim([0 250])
legend('可转移负荷','实际用电计划');
xlabel('时间/h')
ylabel('功率/kw')

figure(6)
bar(result_PG,0.7,'b')
xlim([1 2400])
ylim([0 5000])
% axis([1,24 0 1000])
legend('燃气轮机出力');
xlabel('时间/h')
ylabel('功率/kw')


figure(7)
grid
plot(Cday,'r','linewidth',1)
xlim([1 100])
ylim([-2000 10000])
legend('两节鲁棒优化微网运行总成本');
xlabel('时间/h')
ylabel('总成本/元（RMB）')