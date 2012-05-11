%% Loading dataset.
clear
clc
% Adding plotting paths for myaa
addpath('C:\Users\zigg\Dropbox\Matlab\myaa')
addpath('M:\Dropbox\Matlab\myaa')
addpath('/mn/felt/u9/sigve/Dropbox/Matlab/myaa')
%%
% Importing data
imported_data = importdata('VMC_SGD.dat'); 

thermal = 200;
alpha = imported_data(:,1);
beta = imported_data(:,2);

subplot(2,1,1); 
plot(alpha);
xlabel('\alpha');
hold on

cumAlpha = zeros(length(alpha),1);
for i=1:length(alpha)
   cumAlpha(i) = sum(alpha(thermal:i))/length(alpha(thermal:i));
end

plot(cumAlpha, 'color', [rand rand rand])
hold off
   
subplot(2,1,2);
plot(beta);
xlabel('\beta');

hold on
cumBeta= zeros(length(beta),1);
for i=1:length(beta)
   cumBeta(i) = sum(beta(thermal:i))/length(beta(thermal:i));
end

plot(cumBeta, 'color', [rand rand rand])
hold off
   
a = sum(alpha(thermal:length(alpha)))/length(thermal:length(alpha))
b = sum(beta(thermal:length(beta)))/length(thermal:length(beta))