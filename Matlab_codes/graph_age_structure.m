% plot data from age matrices
ageclass = {'0-4','5-9','10-14','15-19','20-24','25-29','30-34','35-39', ...
    '40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79','80-84',...
    '84-89','90-94','95-99','100+'};
X = categorical(ageclass);
X = reordercats(X,ageclass);
figure
h = barh(X,[F M],3.,'histc');
set(h, {'DisplayName'}, {'Female','Male'}')
set(gca,'Fontsize',20)
xlabel('Age')
