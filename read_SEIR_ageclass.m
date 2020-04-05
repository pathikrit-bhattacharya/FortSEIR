% Read data form unformatted binary file
fname = [pwd '\out\pred_out'];
N_tot = 1.37e9;
fid   = fopen(fname);
nx    = 16*6;
S     = [];
E     = [];
Is     = [];
Ia     = [];
R     = [];
N     = [];
nt    = 3650;
dt    = 0.1;
t     = cumsum(dt*ones(3650,1));
t     = [0;t(1:end-1)];
for i=1:nt
    dat = fread(fid,nx,'double');
    dat = dat';
    S   = [S;dat(1:6:nx)];
    E   = [E;dat(2:6:nx)];
    Is  = [Is;dat(3:6:nx)];
    Ia  = [Ia;dat(4:6:nx)];
    R   = [R;dat(5:6:nx)];
    N   = [N;dat(6:6:nx)];
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
bar_stack = [E(plot_ind,:)' Is(plot_ind,:)' Ia(plot_ind,:)' R(plot_ind,:)'];
ageclass = {'0-4','5-9','10-14','15-19','20-24','25-29','30-34','35-39', ...
    '40-44','45-49','50-54','55-59','60-64','65-69','70-74','75+'};
X = categorical(ageclass);
X = reordercats(X,ageclass);
h = bar(X,log10(bar_stack),1);
set(h, {'DisplayName'}, {'E','I_s','I_a','R'}')
set(gca,'Fontsize',20)
xlabel('Age')
ylabel('log # of people')
fclose(fid);