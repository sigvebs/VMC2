%% Loading dataset.
clear
clc
%% Adding plotting paths for myaa
addpath('C:\Users\zigg\Dropbox\Matlab\myaa')
addpath('M:\Dropbox\Matlab\myaa')
addpath('/mn/felt/u9/sigve/Dropbox/Matlab/myaa')

%%
imported_data = importdata('density.dat'); 

% Storing imported data in the right structures.
N = imported_data(1,1);
dr = imported_data(1,2);
MC_samples = imported_data(1,3);
rho = zeros(N,N);

n = 2;
for i=1:N
   for j=1:N
      rho(i,j) = imported_data(n, 3);
      n = n + 1;
   end;
end;

[x,y] = meshgrid(-N/2*dr:dr:N/2*dr-dr);

%% New plotting style
clc
fig =  figure('Color',[1 1 1]);
logoax = axes('Visible','on', 'parent', fig,  'fontsize', 12);
pl = surf(x, y, rho, 'parent', logoax);

%% Plotting contour.
clc
fig =  figure('Color',[1 1 1]);
logoax = axes('Visible','on', 'parent', fig,  'fontsize', 12);
pl = contourf(x, y, rho,'DisplayName','\rho');
%figure(gcf)
axis equal

%%
colormap('summer')
%colorbar
brighten(0.5)

hYLabel = ylabel('y');
hXLabel = xlabel('x');
hZLabel = zlabel('One Body Density');

grid off

set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde');
set([gca]             , ...
    'FontSize'   , 10           );
set([hXLabel, hYLabel, hZLabel]  , ...
    'FontSize'   , 12           , ...
    'FontWeight' , 'bold'      );

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...
  'ZColor'      , [.1 .1 .1], ...
  'LineWidth'   , 1.2         );

set(gcf, 'PaperPosition', [0 0 12 10]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [12 10]); %Set the paper to have width 5 and height 5.
saveas(gcf, '2particles_n4_w001_DENSITY', 'pdf') %Save figure
%% Smoothing the plot
for I = 1:N
    rho(:,I) = smooth(rho(:,I),'sgolay',2);
end

%%
myaa('publish');

