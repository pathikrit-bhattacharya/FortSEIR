% Code to generate age-classified parameters for the SEIR model
% The parameters are: 
% rho   = probability that an exposed person turns infectious
% kappa = inverse of the mean incubation period (1/t_incubation)
% gamma = inverse of the mean infectious period (1/t_infectious)
% mu_n  = age-specific mortality rate under normal circumstances
% mu_d  = age-specific mortality rate under disese
% w_s = factor controlling effect of lockdown in schools
% w_w = factor controlling effect of lockdown at work
% w_h = factor controlling effect of lockdown at home
% w_o = factor controlling effect of lockdown at other locations
Incub_dur = 6.4;
Infec_dur = 7;
lambda  = 0.018/365; % Rate of transmission of disease
t_tot   = 365; % How long will the simulation run
N_tot   = 1.e4; % Total poulation
rho = [0;0.4*ones(5,1);0.8*ones(10,1)];
alpha = 0.25; % Prob that infected asymptomatic infects
beta  = 0.5; % Prob. that infectious transmits disease to susceptible
kappa = ones(16,1)/Incub_dur;
gamma = ones(16,1)/Incub_dur;
mu_n  = zeros(16,1); % 0.007*ones(16,1)/365; % scale to per day
mu_d  = zeros(16,1); % 0.04*[ones(8,1);2*ones(2,1);6*ones(2,1);18*ones(2,1);57*ones(2,1)]/365; % From worldometer proportions of % death by age (March 31st, 2020)
w_s = 0*ones(16,1); % [1.5*ones(4,1);1.1*ones(12,1)]; % Effect of lockdown at school
w_w = [ones(4,1);0.5*ones(10,1);ones(2,1)]; % Effect of lockdown at work
w_h = [1.5*ones(4,1);1.1*ones(10,1);ones(2,1)]; % Effect of lockdown at home
w_o = [0.5*ones(14,1);ones(2,1)]; % Effect of lockdown at other places
Nis  = load([pwd '\input\India_age_data.txt']); % Population fractions in age groups
I0   = 200*Nis; % Number of infected initially, populated as population age fractions
S0  = Nis*N_tot-200; % Initial number of susceptible in age group
A = [rho,kappa,gamma,mu_n,mu_d,w_s,w_w,w_h,w_o,Nis,I0,S0];
t_l1    = 7; % Lockdown onset
t_l2    = 28; % Lockdown stop
A_last = zeros(16,1); % Last column of A
A_last(1:7) = [alpha;beta;lambda;t_tot;N_tot;t_l1;t_l2]; % The first four values of A_last
A = [A A_last];
fname = [pwd '\input\inparams_SEIR_ageclass'];
fid = fopen(fname,'w+');
fprintf(fid,[repmat('%16.9e ',1,12) '%16.9e\n'],A');
fclose(fid);

