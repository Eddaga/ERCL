%% explain section
% This Code based on "Power Consumption Characterization, Modeling and Estimation of Electric Vehicles"
% energy = integral(IVdt) = integral(Fds)
% I , V -> EV motor Voltage and Current.
% F     -> traction force, s-> moving distance

% P = F * (ds/dt) = F*v = F_{V} = F_{R} + F_{A} + F_{G} + F_{I} + F_{B} 

%where
%F_{R} : rolling     resistance                     ∝ C_{rr}*W
%F_{A} : aerodynamic resistance                     ∝ 1/2 *ρC_{d}Av2
%F_{G} : gradient    resistance                     ∝ W sinθ
%F_{I} : inertia     resistance                     ∝ ma
%F_{B} : brake force proveded by hydraulic brake

% C_{rr} : rolling resistance coefficient
% W      : vehicle weight
% ρC_{d} : air density
% A      : drag coeeicient
% v      : vehicle speed
% θ      : gradient angle
% m      : vehicle mass
% a      : vehicle acceleration

% in reference thesis, Fa is ignored because of vehicle's limited speed(34km/h)
% but our vehicle is for high speed racing, so Fa must considered.

% So, F_{V} ≈ (α+βsinθ+γa+)mv + δv^2 

% α is rolling     resistance
% β is gradient    resistance
% γ is inertia     resistance
% δ is aerodynamic resistance

% (α+βsinθ+γa)v + δv^2

% θ -> in frontyard drive, it's always 0
% m is coefficient
% so, the P equation is
%% 
% P = (α+γa)v + δv^2

% to LR, v as v_! make av as x_2, v^2 as x_3 => equation is 
% P = α*x_1 + γ*x_2 + δ*x_3
% => all degrees are linear.

% velocity is motor's velocity, so we have to LR for velocity(v) first.

% SO, work like this.

%%
%1. load Data
%    1.1 load Data V_1, V_2, V_3, C
%    1.2 laod motor velocity
%2. calculate loaded data for P, a, v
%    2.1 get P. P = C * (V_1+V_2+V_3)
%    2.2 get a. a = calculas of v
%3.  make x_1, x_2, x_3


%% working section
%% 1. load Data
    data = readtable("data.xlsx");
%%
    %1.1 load Data V_1, V_2, V_3, C
    tlength = height(data(:,"frontV"));
    V1 = table2array(data(2500:tlength,"frontV"));
    V2 = table2array(data(2500:tlength,"backV"));
    V3 = table2array(data(2500:tlength,"trunkV"));
    C = table2array(data(2500:tlength,"frontC"));

%   1.2 laod motor velocity
%%
    motorVelocity = table2array(data(2500:tlength,"aMotorVelocity"));

%% 2. calculate loaded data for P, a, v
    %2.1 get P. P = C * (V_1+V_2+V_3)
    P = (V1 + V2 + V3).*C;
    
    %2.2 get a. a = calculas of v
    dmv = diff(motorVelocity);
    a = dmv./0.01;
    
%% 3. make x_1, x_2, x_3 and do linear Regression 
    P = P(1:length(P)-1,1);
    motorVelocity = motorVelocity(1:length(motorVelocity)-1,1);
    variables = [motorVelocity, a.*motorVelocity, motorVelocity.^2];% ones(size(motorVelocity))];
    coefficients = variables \ P;
    
%% 4. compare with machin learing result.
test_target = readmatrix('dataframes.xlsx','Sheet','test_target');
test_output = readmatrix('dataframes.xlsx','Sheet','test_outputs');

test_data = [test_target(:,2) test_output(:,2)];

for i = 1 : 1 : length(test_data)
    test_data(i,3) = (test_data(i,1) - test_data(i,2)) / test_data(i,1);
end

xAxis = 1: 1: length(test_data)-1;

%% make LRd data
test_input = readmatrix('test_data.csv');
test_motorVelocity = test_input(:,2);

test_dmv = diff(test_motorVelocity);
test_a = test_dmv./0.01;

lrdP = (coefficients(1) .*test_motorVelocity(1:length(test_a))) + ...
       (coefficients(2) *test_a .* test_motorVelocity(1:length(test_a))) + ...
       (coefficients(3) .* (test_motorVelocity(1:length(test_a)).^2)) ; % + ...
       %coefficients(4);

%% plot
%figure(1)
%plot(xAxis,test_data(:,3));

mvdTarget = movmean(test_data(:,1),100);
mvdOutput = movmean(test_data(:,2),100);
mvdLrdOutput = movmean(lrdP,100);
figure(3)
plot(mvdLrdOutput);
figure(2)
mvAveragedGraph = nexttile;
plot(mvAveragedGraph,xAxis,mvdTarget(1:27981,1),xAxis,mvdOutput(1:27981,1));
title(mvAveragedGraph,'LearnResult');
xlabel(mvAveragedGraph,'tick');
ylabel(mvAveragedGraph,'P(V*I)');
legend(mvAveragedGraph,'target','output');

