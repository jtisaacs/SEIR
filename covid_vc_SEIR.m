clearvars;close all;clc;
[tableConfirmed,tableDeaths,tableRecovered,time] = getDataCOVID_US();
timeRef = time;

Location = 'Ventura, California, US'; % Find every cities in Washington state
% Location = 'New York'; % Find every cities in New York state
try
    indC = find(contains(tableConfirmed.Combined_Key,Location)==1);
    indD = find(contains(tableDeaths.Combined_Key,Location)==1);
catch exception    
    searchLoc = strfind(tableConfirmed.Combined_Key,Location);
    indC = find(~cellfun(@isempty,searchLoc))  ;
    
    searchLoc = strfind(tableDeaths.Combined_Key,Location);
    indD = find(~cellfun(@isempty,searchLoc))   ; 
end

disp(tableConfirmed(indC,11));

% disp(tableConfirmed(indC,1:2))


% Initialisation
Confirmed = 0;
Deaths = 0;
Npop = 0;
for ii=1:numel(indC)
    Confirmed =  Confirmed + table2array(tableConfirmed(indC(ii),12:end));
end
for ii=1:numel(indD)
    Deaths =  Deaths + table2array(tableDeaths(indD(ii),13:end));
    Npop= Npop + table2array(tableDeaths(indD(ii),12)); % population (dummy number here)
end

%%%% JTI is adding in data from VC Emergencey for Confirmed, Deaths, and
%%%% Recovered

%% Data From Erin
Recovered = [0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
1
1
2
3
4
6
6
8
8
10
20
24
35
46
57
59
62
87
101
113
124
147
150
152
175
191
202
219
240
242
242
260
280
286
306
327
335
358
375
383
390
401
402]';

Deaths = [0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
1
1
1
2
2
4
4
5
5
5
7
8
8
8
8
8
9
9
12
13
14
14
14
14
14
14
14
14
14
16
16
16
16
16
17
17
17
17
19
19
19]';

Confirmed = [1
1
1
1
1
3
3
4
8
11
13
15
25
34
48
61
70
74
81
103
125
147
157
173
178
187
206
216
234
246
265
274
278
297
315
329
346
362
368
374
391
408
425
437
451
455
459
474
487
499
508
515
518
521
542
554
563
574
577]';

t0 = datetime(2020,3,1); % The date of the first confirmed case was March 1, 2020.
timeRef = datetime(t0:1:datetime(datestr(floor(datenum(t0))+datenum(size(Confirmed,2)-1))));
days = 0:1:length(time)-1;


% If the number of confirmed Confirmed cases is small, it is difficult to know whether
% the quarantine has been rigorously applied or not. In addition, this
% suggests that the number of infectious is much larger than the number of
% confirmed cases
time = timeRef;
minNum= round(0.1*max(Confirmed)); % 5% of the maximal number of confirmed is used for the initial conditions
Recovered(Confirmed<=minNum)=[];
Deaths(Confirmed<=minNum)=[];
time(Confirmed<=minNum)= [];
time(1)
days(Confirmed<=minNum)= [];
Confirmed(Confirmed<=minNum)=[];

fprintf(['Population = ',num2str(Npop),' \n'])

% Definition of the first estimates for the parameters
alpha_guess = 0.05;
beta_guess = 0.8; % Infection rate
LT_guess = 5; % latent time in days
Q_guess = 0.5; % rate at which infectious people enter in quarantine
lambda_guess = [0.1,0.1,10]; % recovery rate
kappa_guess = [0.01,0.01,10]; % death rate

guess = [alpha_guess,...
    beta_guess,...
    1/LT_guess,...
    Q_guess,...
    lambda_guess,...
    kappa_guess];

E0 = Confirmed(1); % Initial number of exposed cases. Unknown but unlikely to be zero.
I0 = Confirmed(1); % Initial number of infectious cases. Unknown but unlikely to be zero.
Q0 = Confirmed(1)-Deaths(1)-Recovered(1);
R0 = Recovered(1); 
D0 = Deaths(1);

% Parameter estimation with the lsqcurvefit function[alpha1,beta1,gamma1,delta1,Lambda1,Kappa1] = ...
    [alpha_s,beta_s,gamma_s,delta_s,Lambda_s,Kappa_s,lambdaFun_s,kappaFun_s] = ...
        fit_SEIQRDP(Confirmed(1:17)-Deaths(1:17)-Recovered(1:17),Recovered(1:17),Deaths(1:17),Npop,E0,I0,time(1:17),guess,'Display','off');

% Parameter estimation with the lsqcurvefit function[alpha1,beta1,gamma1,delta1,Lambda1,Kappa1] = ...
    [alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,lambdaFun,kappaFun] = ...
        fit_SEIQRDP(Confirmed-Deaths-Recovered,Recovered,Deaths,Npop,E0,I0,time,guess,'Display','off');
    
    
    
dt = 1/24; % time step

% time1 = datetime(time(1)):dt:datetime(datestr(floor(datenum(now))));
% time_proj = datetime(time(1)):dt:datetime(datestr(floor(datenum(now))+datenum(10)));
% N = numel(time1);
% t = [0:N-1].*dt;

time1 = time(1):dt:time(end);%datetime(time(1)):dt:datetime(time(end));
N = numel(time1);
t = [datenum(time(1))-datenum(t0):datenum(time(1))-datenum(t0)+N-1].*dt;

time_proj = time(1):dt:datetime(datestr(time(end)+datenum(60)));
N2 = numel(time_proj);
t_proj = [datenum(time(1))-datenum(t0):datenum(time(1))-datenum(t0)+N2-1].*dt;

time_u = time_proj(length(time1):end);
t_u = t_proj(length(t):end);

[S,E,I,Q,R,D,P] = SEIQRDP(alpha1,beta1,...
    gamma1,delta1,Lambda1,Kappa1,Npop,[],E0,I0,Q0,R0,D0,[],t,lambdaFun,kappaFun);

lambda = lambdaFun(Lambda1,t);
% R0_orig = beta1./lambda(1)

[Sp,Ep,Ip,Qp,Rp,Dp,Pp] = SEIQRDP(alpha1,beta1,...
    gamma1,delta1,Lambda1,Kappa1,Npop,[],E0,I0,Q0,R0,D0,[],t_proj,lambdaFun,kappaFun);

lambda = lambdaFun(Lambda1,t_proj);
% R0_proj = beta1./lambda(1)

beta_s; % This is the beta observed in VC up to March 17 when StayAtHome order initiated.
% beta_initial = 0.5; %Lets say we only go back to half interactions.
[Su,Eu,Iu,Qu,Ru,Du,Pu] = SEIQRDP(alpha1,beta_s,...
    gamma1,delta1,Lambda1,Kappa1,Npop,S(end),E(end),I(end),Q(end),R(end),D(end),P(end),t_u,lambdaFun,kappaFun);

lambda = lambdaFun(Lambda1,t_u);
% R0_u = beta1./lambda(1)


checkRates(time,Confirmed-Deaths-Recovered,Recovered,Deaths,kappaFun,lambdaFun,Kappa1,Lambda1);
%%
figure
semilogy(time_proj,Qp+Rp+Dp,'b',time_proj,Rp,'r',time_proj,Qp,'g',time_proj,Dp,'k');
hold on
semilogy(time_u,Qu+Ru+Du,'b',time_u,Ru,'r',time_u,Qu,'g',time_u,Du,'k');
semilogy(time,Confirmed,'bo',time,Recovered,'ro',time,Confirmed-Recovered-Deaths,'go',time,Deaths,'ko');
% ylim([0,1.1*Npop])
ylabel('Number of cases')
xlabel('time (days)')
leg = {'Tested positive  (fitted)','Recovered (fitted)',...
    'Active(fitted)','Deceased (fitted)',...
    'Tested positive (reported)',...
    'Recovered (reported)','Active (reported)','Deceased  (reported)'};
legend(leg{:},'location','southoutside')
set(gcf,'color','w')
grid on
axis tight
title([Location])
set(gca,'yscale','lin')