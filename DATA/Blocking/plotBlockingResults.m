%% Loading dataset.
clear
clc
addpath('C:\Users\zigg\Dropbox\Matlab\myaa')
addpath('M:\Dropbox\Matlab\myaa')

imported_data = importdata('blockingResults.dat');
% Storing imported data in the right structures.
N = length(imported_data(:,1));
dN = imported_data(2,1) - imported_data(1,1);
E = imported_data(:,2);
sigma = imported_data(:,3);
blocks = imported_data(:,1);
exact_E = zeros(N,1) + 3.0;

myMap = autumn(3);
% Plotting error
clc
color = [177, 0, 38]/256;
fig =  figure('Color',[1 1 1]);
pl = plot(blocks, sigma,'.','color', color, 'LineWidth',0.5);%, 'linesmoothing', 'on', 'color', myMap(2,:), 'LineWidth',1.5);

hYLabel = ylabel('\sigma');
hXLabel = xlabel('Block size');

grid off

set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde');
set([gca]             , ...
    'FontSize'   , 10           );
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 12          , ...
    'FontWeight' , 'bold'      );

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         );

%%
set(gcf, 'PaperPosition', [0 0 12 10]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [12 10]); %Set the paper to have width 5 and height 5.
saveas(gcf, '2particles_n6_w1_BF', 'pdf') %Save figure
%%
myaa('publish');
