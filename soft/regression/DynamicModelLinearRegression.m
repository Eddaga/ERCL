%% explain section
% This Code based on "Power Consumption Characterization, Modeling and Estimation of Electric Vehicles"
% energy = integral(IVdt) = integral(Fds)
% I , V -> EV motor Voltage and Current.
% F     -> traction force, s-> moving distance

% P = F * (ds/dt) = F_{V} = F_{R} + F_{A} + F_{G} + F_{I} + F_{B} 

%where
%F_{R} : rolling     resistance                     ∝ C_{rr}*W
%F_{A} : aerodynamic resistance                     ∝ 1/2 *ρC_{d}Av2
%F_{G} : gradient    resistance                     ∝ W sinθ
%F_{I} : inertia     resistance                     ∝ ma
%F_{B} : brake force proveded by hydraulic brake

% C_{rr} : rolling resistance coefficient
% W      : vehicle weight
% C_{d}  : air density
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
% so, the P equation is
% P = (α+γa)v + δv^2

%% working section
    %% data load
data = readtable("data.xlsx");
tlength = height(data(:,"frontV"));
V1 = table2array(data(2500:tlength,"frontV"));
V2 = table2array(data(2500:tlength,"backV"));
V3 = table2array(data(2500:tlength,"trunkV"));
C = table2array(data(2500:tlength,"frontC"));
P = (V1 + V2 + V3).*C;


aXx = table2array(data(2500:tlength,"AccX"));
aYy = table2array(data(2500:tlength,"AccY"));
aZz = table2array(data(2500:tlength,"AccZ"));
gXx = table2array(data(:,"GyroX"));
gYy = table2array(data(:,"GyroY"));
gZz = table2array(data(:,"GyroZ"));

v = table2array(data(2500:tlength,"aMotorVelocity"));
%% get accel
% Define constants
g = 9.81; % gravitational acceleration

aX = aXx * g;
aY = aYy * g;
aZ = aZz * g;
gX = gXx * pi / 180.0;
gY = gYy * pi / 180.0;
gZ = gZz * pi / 180.0;
%{
% Calculate pitch and roll angles
pitch = atan2(aX, sqrt(aY.^2 + aZ.^2));
roll = atan2(aY, sqrt(aX.^2 + aZ.^2));

% Calculate actual angular velocity
gX_actual = gX - gY.*sin(roll) + gZ.*cos(roll);
gY_actual = gY.*cos(pitch) + gZ.*sin(pitch);
gZ_actual = -gY.*sin(pitch) + gZ.*cos(pitch);

% Integrate actual angular velocity over time to get change in orientation
dt = 0.01; % sample time (in seconds)
delta_pitch = cumsum(gX_actual * dt);
delta_roll = cumsum(gY_actual * dt);
delta_yaw = cumsum(gZ_actual * dt);

% Apply rotation matrix to accelerometer readings to get actual linear acceleration
accX_actual = aX.*cos(delta_pitch) + aY.*sin(delta_roll).*sin(delta_pitch) + aZ.*cos(delta_roll).*sin(delta_pitch);
accY_actual = aY.*cos(delta_roll) - aZ.*sin(delta_roll);
accZ_actual = -aX.*sin(delta_pitch) + aY.*sin(delta_roll).*cos(delta_pitch) + aZ.*cos(delta_roll).*cos(delta_pitch);

% Combine linear acceleration and gravitational acceleration to get total acceleration
a = sqrt(accX_actual.^2 + accY_actual.^2 + accZ_actual.^2) - g;
%}
dt = 0.01;
a = sqrt(aX.^2+aY.^2+aZ.^2)-g ;
minaZz = min(aZz);
minaZ = min(aZ);
movemeanda = movmean(a(:,1),100);
movemeandat = linspace(0, length(movemeanda)*dt,length(movemeanda));
plot(movemeandat,movemeanda);
figure;
xlabel('Time (s)');
ylabel('AccelerationMovemeand (m/s^2)');
% Plot acceleration data
t = linspace(0, length(a)*dt, length(a));
plot(t, a);
figure;
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
%%

%% perform linear regression
X = [a.*v, v.^2, ones(length(a),1)];
lm = fitlm(X, P);
coef = lm.Coefficients.Estimate;
ci = lm.coefCI;

%% plot results
figure;
scatter3(a,v,P, 'filled');
hold on;
a_lim = [min(a), max(a)];
v_lim = [min(v), max(v)];
[X1,X2] = meshgrid(a_lim,v_lim);
Y = coef(1)*X1.*X2 + coef(2)*X2.^2 + coef(3);
surf(X1,X2,Y,'FaceAlpha',0.5);
xlabel('Acceleration');
ylabel('Velocity');
zlabel('Power');
title('Linear Regression Model');