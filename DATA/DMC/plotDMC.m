%% Loading dataset.
clear
clc
imported_data = importdata('dmc.dat'); 

thermal = 20;
E = imported_data(:,1);

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

Etot = sum(E(thermal:length(E)))/length(thermal:length(E))