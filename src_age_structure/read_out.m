% File to read output of the SEIR epidemic simulation
% First line: Header --> alpha, beta, gamma, rho, R0 --> First line of values are model parameters
% Second line: Format >>> t S(t) E(t) I(t) R(t)
nt      = 36500;
N_tot   = 1e9;
outfile = [pwd '\out\pred_out'];
fid     = fopen(outfile,'r+');
params  = textscan(fid,'%f %f %f %f %f',1);
dat     = textscan(fid,'%f %f %f %f %f',nt);
t_s     = dat{1};
s_t     = dat{2};
e_t     = dat{3};
i_t     = dat{4};
r_t     = dat{5};
figure
plot(t_s,e_t,'linewidth',2,'DisplayName','Exposed')
hold on
plot(t_s,i_t,'linewidth',2,'DisplayName','Infected')
legend
xlabel('Time [days]')
ylabel('# of People')
set(gca,'Fontsize',18)
title(['Model params: t_{inc}^{-1}=' num2str(params{1},'%5.3f') ', \beta=' num2str(params{2},'%5.3f'),...
    ', R_0=' num2str(params{5},'%5.3f'), ', t_{d}^{-1}=' num2str(params{3},'%5.3f')])
fclose(fid);