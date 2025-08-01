function [obj1, x, y] = Primal_Problem(k,uncert,price,PDR_exp)
%Initialize fixed constants
PGmax = 800;
PGmin = 0;
a = 0.67;
b = 0;
PSmax = 500;
ESmax = 1800;
ESmin = 400;
ES_0 = 1000;
Ks = 0.38;
eta = 0.95;
KDR = 0.32;
DDR = 2940;
dDR = 0.4;
PMmax = 1000;
gammaPV = 6;
gammaL = 12;



% %Initializing the grid model
% mpc = loadcase('case39_ED.m');
% mpc_load = load('load3996.mat');
% 
% load = mpc_load.Pd;
% bus_num = length(mpc.bus(:,1));
% gen_num = length(mpc.gen(:,1));
% storage = ex_storage(mpc);
% wind = ex_wind_uc(mpc);
% [i_ess,mpc] = addstorage(storage,mpc);
% [i_wind,mpc] = addwind(wind,mpc);
% 
% Pbranch_max = mpc.branch(:,6);
% Pbranch_min = -Pbranch_max;
% 
% HP = makePTDF(mpc);

% %Electricity price per hour
% price =[0.47,0.47,0.47,0.47,0.47,0.47,0.47,0.9,1.35,1.35,1.35,0.9,0.9,0.9,0.9,0.9,0.9,0.9,1.35,1.35,1.35,1.35,1.35,0.47];
% 
% %Estimated variable load 
% PDR_exp=[80,70,60,50,70,70,90,100,120,150,170,200,140,100,100,120,140,150,190,200,200,190,100,80];
% PDR_exp_100D = repmat(PDR_exp, 1, 100);

M=1e6;


n_time = 24;
delta_t = 1;



%Decide variables
bus_p = sdpvar(1,n_time,'full');
PG = sdpvar(1,n_time,'full');
PSchg  = sdpvar(1,n_time,'full');
PSdis  = sdpvar(1,n_time,'full');
US  = sdpvar(1,n_time,'full');
PDR = sdpvar(1,n_time,'full');
PDR1 = sdpvar(1,n_time,'full');
PDR2 = sdpvar(1,n_time,'full');
DDRmin = sdpvar(1,n_time,'full');
DDRmax = sdpvar(1,n_time,'full');
Psell = sdpvar(1,n_time,'full');
Pbuy = sdpvar(1,n_time,'full');
UM = sdpvar(1,n_time,'full');
PPV=sdpvar(1,n_time,'full');
PL=sdpvar(1,n_time,'full');

%Dual problem's variables
gamma=sdpvar(12*n_time,1,'full');
lambda=sdpvar(50,1,'full');
miu=sdpvar(4*n_time,1,'full');
phi=sdpvar(2*n_time,1,'full');

% Boolean/binary variables
Bpv=binvar(1,n_time,'full');
BL=binvar(1,n_time,'full');
B=[Bpv,BL]';
B1=sdpvar(2*n_time,1,'full');



%Setting up coefficient and constant matrix

%Objective function and double coefficient matrix
x=[US,UM]';
y=[PG,PSchg,PSdis,PDR,PDR1,PDR2,Pbuy,Psell,PPV,PL]';
c=[a*ones(1,24),Ks*eta*ones(1,24),Ks/eta*ones(1,24),zeros(1,24),KDR*ones(1,24),KDR*ones(1,24),price,-price,zeros(1,48)];

%Inequality constraints coefficient and constant matrix
D=[-eye(24,240);...
    eye(24,240);...
    zeros(24,24),-eta*tril(ones(24,24),0),(1/eta)*tril(ones(24,24),0),zeros(24,168);...
    zeros(24,24),eta*tril(ones(24,24),0),-(1/eta)*tril(ones(24,24),0),zeros(24,168);...
    zeros(24,72),-eye(24,24),zeros(24,144);...
    zeros(24,72),eye(24,24),zeros(24,144);...
    zeros(24,24),eye(24,216);...
    zeros(24,48),eye(24,192);...
    zeros(24,96),eye(24,144);...
    zeros(24,120),eye(24,120);...
    zeros(24,144),eye(24,96);...
    zeros(24,168),eye(24,72)];
d=[-PGmax*ones(1,24),PGmin*ones(1,24),-(ESmax-ES_0)*ones(1,24),(ESmin-ES_0)*ones(1,24),-PDR_exp*(1+dDR),PDR_exp*(1-dDR),zeros(1,144)]';

%Equality constraints coefficient and constant matrix
K=[zeros(1,24),eta*ones(1,24),(-1/eta)*ones(1,24),zeros(1,168);...
    zeros(1,72),ones(1,24),zeros(1,144);...
    zeros(24,72),eye(24,24),eye(24,24),-eye(24,24),zeros(24,96);...
    eye(24)     -eye(24)    eye(24)     -eye(24)    zeros(24,48)    eye(24) -eye(24)    eye(24) -eye(24)];

%Vector matrix
g=[0,DDR,PDR_exp,zeros(1,24)]';

%Double variable inequality constraints coefficient and constant matrixes
F=[-PSmax*eye(24,48);PSmax*eye(24,48);zeros(24,24),PMmax*eye(24,24);zeros(24,24),-PMmax*eye(24,24)];
G=[zeros(24,24),-eye(24,216);zeros(24,48),-eye(24,192);zeros(24,144),-eye(24,96);zeros(24,168),-eye(24,72)];
h=[-PSmax*ones(1,24),zeros(1,24),zeros(1,24),-PMmax*ones(1,24)]';


% Undetermined values/fluctuations coefficient matrix
I=[zeros(24,192),eye(24,48);...
    zeros(24,216),eye(24,24) ];
% u=[Ppv,PL]';


%Setting up constraints 
cons=[];
u=sdpvar(48,1,'full');
%Primal Problem 
cons=[cons,D*y>=d];
cons=[cons,K*y==g];
cons=[cons,F*x+G*y>=h];
cons=[cons,I*y==uncert];

% %Dual Problem
% cons=[cons,D'*gamma+K'*lambda+G'*miu+I'*phi<=c'];
% cons=[cons,gamma>=0,miu>=0];
% cons=[cons,u==u_f+du_f.*B];
% cons=[cons,sum(Bpv)<=gammaPV];
% cons=[cons,sum(BL)<=gammaL];
% 
% 
% cons=[cons,-M*B<=B1,B1<=M*B];
% cons=[cons,du_f.*phi-M*(1-B)<=B1,B1<=du_f.*phi];



%Setting the Problem
obj1 = c*y;
ops = sdpsettings('solver','gurobi','verbose',0);
res = optimize(cons,obj1,ops); 
if res.problem == 0

else
    fprintf('The %d th iteration is unsuccessfull\n',k);
end


obj1 = value(obj1);
x = value(x);
y = value(y);

result_y = double(y);
PG_res=result_y(1:24,1);
PSchg_res=result_y(25:48,1);
PSdis_res=result_y(49:72,1);
PDR_res=result_y(73:96,1);
Pbuy_res=result_y(145:168,1);
Psell_res=result_y(169:192,1);


end
