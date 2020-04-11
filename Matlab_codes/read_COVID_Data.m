% Extract state level or India level data from covid_data_Apr8.mat
clear
load([pwd '\input\covid_data_Apr8.mat'])
state_list = unique(StateUnionTerritory); % States which have had reported cases
dat = struct();
dat.state = state_list;
dat.Date  = unique(Date);
dat.CumCount = zeros(length(dat.Date),length(state_list));
date.Population = zeros(32,1);
for i = 1:length(Date)
    r = find(StateUnionTerritory(i)==dat.state);
    q = find(Date(i)==dat.Date);
    dat.CumCount(q,r) = Confirmed(i);
end
load([pwd '\input\India_census_2011.mat'])
for i = 1:length(pop.state)
    q = find(pop.state(i)==dat.state);
    dat.Population(q) = pop.Populaion(i);
end
state = categorical(cellstr('Odisha'));
n = find(dat.state==state);
%plot(dat.Date,dat.CumCount(:,n))
date = dat.Date; I_India = sum(dat.CumCount,2);day=1:length(dat.Date);
plot(date,I_India)
save ./input/India_ts day date I_India