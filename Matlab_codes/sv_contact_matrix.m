% File to read out country specific contact matrix data from
% Projecting social contact matrices in 152 countries using contact surveys and
% demographic data, Kiesha Prem, A. R. Cook, Mark Jit, PLOS Comp. Biol.,
% 2017
datadir = [pwd '\contact_matrices_152_countries'];
dir_content = dir(datadir);
sz = length(dir_content);
k =  0;
cmap = jet(20);
cmap = flipud(cmap(1:10,:));
xp = [0.05 0.5];
hp = [0.025 0.375 0.725];
w  = 0.4;
h  = 0.25;
pos = [];
% Make grid for subplots
for i = 1:3
    for j = 1:2
        pos = [pos;xp(j) hp(i) w h]; %  [left bottom width height]
    end
end
figure
% Actual plots
for i = 3:2:12
   name = [datadir '\' dir_content(i).name];
   cm_type = dir_content(i).name;
   A = xlsread(name,'Germany','A2:P17');
   c = i-2-k;
   k = k+1; 
   subplot('Position',pos(k,:))
   imshow(A,'Colormap',cmap,'XData',1:80,'YData',1:80)
   lab_plot = cm_type(13:end-5);
   lab_plot(lab_plot=='_'|lab_plot=='1')=' ';
   ylabel(lab_plot)
   axis on; box on
   set(gca,'XAxisLocation','top')
   xticks(0:10:80)
   yticks(0:10:80)
   lab_plot(lab_plot==' ')='_';
%    fname = [pwd '\input\mat_files_India\' lab_plot(1:end-2) '.txt'];
%    fid   = fopen(fname,'w+');
%    fprintf(fid,[repmat('%16.9e ',1,15) '%16.9e\n'],A');
%    fclose(fid);
end
