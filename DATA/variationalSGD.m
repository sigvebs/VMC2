%% Plotting SGD generated parameters and priting average results.
clear
clc

% Importing data
imported_data = importdata('VMC_SGD.dat'); 
SGDSamples = imported_data(1,1);
McSamples = imported_data(1,2);
walkers  = imported_data(1,3);
nParticles = imported_data(1,4);
w = imported_data(1,5);

%
thermal = 30;
alpha = imported_data(2:length(imported_data(:,1)),1);
beta = imported_data(2:length(imported_data(:,2)),2);

subplot(2,1,1); 
plot(alpha);
xlabel('\alpha');

subplot(2,1,2);
plot(beta);
xlabel('\beta');

% Calculating average after thermalization
a = sum(alpha(thermal:length(alpha)))/length(thermal:length(alpha));
b = sum(beta(thermal:length(beta)))/length(thermal:length(beta));

% Caulating variance
varAlpha = dot(alpha(thermal:length(alpha)), alpha(thermal:length(alpha)))/length(thermal:length(alpha)) - a^2;
varBeta = dot(beta(thermal:length(beta)), beta(thermal:length(beta)))/length(thermal:length(beta)) - b^2;

% Writing results in latex format
fprintf('%i & %3d & %4d & %4d & %4d & %4d \\\\', nParticles, w, a, varAlpha, b, varBeta );