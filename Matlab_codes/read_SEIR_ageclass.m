% Read data form unformatted binary file
fname = [pwd '\out\pred_out'];
fRname = [pwd '\out\R_out'];
% what tstep to plot histogram: options are 'max' or 'final'
step_choose = 'max';
N_tot = 1.37e9;
fid   = fopen(fname);
fRid  = fopen(fRname);
nac   = 16;
neq   = 6;
nx    = nac*neq;
S     = [];
E     = [];
Is    = [];
Ia    = [];
R     = [];
N     = [];
R0    = [];
nt    = 365;
dt    = 1.d0;
t     = cumsum(dt*ones(nt,1));
t     = [0;t(1:end-1)];
for i=1:nt
    dat = fread(fid,nx,'double');
    T   = fread(fRid,[3*nac 3*nac],'double');
    Sigma = fread(fRid,[3*nac 3*nac],'double');
    dat = dat';
    S   = [S;dat(1:6:nx)];
    E   = [E;dat(2:6:nx)];
    Is  = [Is;dat(3:6:nx)];
    Ia  = [Ia;dat(4:6:nx)];
    R   = [R;dat(5:6:nx)];
    N   = [N;dat(6:6:nx)];
    lambda = eig(-T*inv(Sigma));
    R0  = [R0;max(lambda)];
end
S_tot = sum(S,2);
Ia_tot = sum(Ia,2);
Is_tot = sum(Is,2);
I_tot =  Ia_tot + Is_tot;
[~,plot_ind] = max(I_tot);
E_tot = sum(E,2);
R_tot = sum(R,2);
Sum_tot = S_tot+E_tot+I_tot+R_tot;
figure
hold on
plot(t,I_tot,'r-','Linewidth',2,'DisplayName','Total Infected')
plot(t,Is_tot,'b-','Linewidth',2,'DisplayName','Infected - symptomatic')
plot(t,Ia_tot,'b--','Linewidth',2,'DisplayName','Infected - asymptomatic')
box on
xlabel('Time from onset [days]')
ylabel('# of People')
set(gca,'Fontsize',20)
figure
if strcmp(step_choose,'max')
    tstep = plot_ind;
elseif strcmp(step_choose,'final')
    tstep = length(I_tot);
end
bar_stack = [E(tstep,:)' Is(tstep,:)' Ia(tstep,:)' R(tstep,:)'];
ageclass = {'0-4','5-9','10-14','15-19','20-24','25-29','30-34','35-39', ...
    '40-44','45-49','50-54','55-59','60-64','65-69','70-74','75+'};
X = categorical(ageclass);
X = reordercats(X,ageclass);
h = bar(X,log10(bar_stack),1);
set(h, {'DisplayName'}, {'E','I_s','I_a','R'}')
set(gca,'Fontsize',20)
xlabel('Age')
ylabel('log # of people')
figure;
h = plot(t,R0);
xlabel('Time from onset [days]')
ylabel('R_{eff}(t)')
set(gca,'Fontsize',20)
fclose(fid);
fclose(fRid);