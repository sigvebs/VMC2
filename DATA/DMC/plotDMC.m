%% Loading dataset.
clear
clc

% Reading data from file
imported_data = importdata('dmc.dat'); 
nParticles = imported_data(1,1);
w = imported_data(1,2);
DMCSamples  = imported_data(1,3);
blockSize = imported_data(1,4);

thermal = 20;
E = imported_data(2:length(imported_data(:,1)),1);

% Ploting block energy
plot(E, 'color', [rand rand rand]);
xlabel('E');
hold on

% Computing cumulative energy average
cumE = zeros(length(E),1);
for i=1:length(E)
   cumE(i) = sum(E(thermal:i))/length(E(thermal:i));
end

% PLotting cumulative energy average
plot(cumE, 'color', [rand rand rand])
hold off

Etot = sum(E(thermal:length(E)))/length(thermal:length(E));

% Variance
var = dot(E(thermal:length(E)), E(thermal:length(E)))/length(thermal:length(E)) - Etot^2;

% Writing results in latex format
fprintf('%i & %.2f & %.6g & %.5g \\\\', nParticles, w, Etot, var);