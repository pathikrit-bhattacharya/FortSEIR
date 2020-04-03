% Read data form unformatted binary file
fname = [pwd '\out\pred_out'];
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
I_tot = sum(Ia,2) + sum(Is,2);
plot(t,I_tot)
fclose(fid);