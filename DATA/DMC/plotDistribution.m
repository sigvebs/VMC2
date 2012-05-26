%% Loading dataset.

% Reading data from file
importedDistribution = importdata('distributionThermal.dat'); 

% DMC thermal
x = importedDistribution(:,1);
y = importedDistribution(:,2);
r = sqrt(x.*y);

% VMC
importedDistributionVMC = importdata('distributionVMC.dat'); 

xVMC = importedDistributionVMC(:,1);
yVMC = importedDistributionVMC(:,2);
rVMC = sqrt(xVMC.*yVMC);

% DMC
importedDistributionDMC = importdata('distributionDMC.dat'); 
xDMC = importedDistributionDMC(:,1);
yDMC = importedDistributionDMC(:,2);
rDMC = sqrt(xDMC.*yDMC);

% PLotting
subplot(3,1,1)
plot(xDMC, yDMC, 'ro');
legend('DMC')
axis equal

subplot(3,1,2)
r_ = 0.0001:0.05:2;

rDMC = histc(rDMC, r_);
plot(rDMC, 'r');
legend('DMC')
hold on

rVMC = histc(rVMC, r_);
plot(rVMC);
legend('DMC', 'VMC')
hold off

subplot(3,1,3)
plot(xVMC, yVMC, '+')
legend('VMC')
axis equal
