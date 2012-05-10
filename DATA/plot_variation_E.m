%% Loading dataset.
clear
clc
% Adding plotting paths for myaa
addpath('C:\Users\zigg\Dropbox\Matlab\myaa')
addpath('M:\Dropbox\Matlab\myaa')
addpath('/mn/felt/u9/sigve/Dropbox/Matlab/myaa')

% Importing data
imported_data = importdata('VMC_prf.dat'); 
%%
% Storing imported data in the right structures.
n_particles = imported_data(1, 1);
alpha_steps = imported_data(1, 2);
beta_steps = imported_data(1, 3);
MC_samples = imported_data(1, 4);
E = zeros(alpha_steps, beta_steps);
E_sq = zeros(alpha_steps, beta_steps);
alpha = zeros(alpha_steps, beta_steps);
beta = zeros(alpha_steps, beta_steps);

n = 2;
for i=1:alpha_steps
   for j=1:beta_steps
      alpha(i,j) = imported_data(n, 1);
      beta(i,j) = imported_data(n, 2);
      E(i,j) = imported_data(n, 3);
      E_sq(i,j) = imported_data(n, 4);
      n = n + 1;
   end;
end;

% Plotting contour.
clc
fig =  figure('Color',[1 1 1]);
logoax = axes('Visible','on', 'parent', fig,  'fontsize', 12);
pl = contourf(alpha, beta, E, 10, 'DisplayName','E');

% Making it pretty
colormap('winter');
brighten(0.6)
c = colorbar;

hXLabel = xlabel('\alpha');
hYLabel = ylabel('\beta');
hZLabel = xlabel(c, 'E');

grid off

set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde');
set([gca]             , ...
    'FontSize'   , 10           );
set([hXLabel, hYLabel, hZLabel]  , ...
    'FontSize'   , 10           , ...
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
saveas(gcf, 'variational_plot', 'pdf') %Save figure

%%
myaa('publish');
