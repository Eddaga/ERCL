%% explain section
% git passwd = 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
% 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
% fucking git cola.. give me linux-cola!!!

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
V1 = data(:,"frontV");
V2 = data(:,"backV");
V3 = data(:,"trunkV");
C = data(:,"frontC");
y = (V1 + V2 + V3) * C;
x1 = data(:,"accel");
x2 = data(:,"aMotorvelocity");

% 데이터 샘플 개수
m = length(y);

% X 변수를 만듭니다.
X = [ones(m, 1), x1, x2, x1.*x2, x2.^2];

% theta(모델 파라미터) 값을 초기화합니다.
initial_theta = zeros(size(X, 2), 1);

% 옵션을 설정합니다.
options = optimset('GradObj', 'on', 'MaxIter', 400);

% 비용 함수와 기울기를 계산하는 함수를 정의합니다.
costFunction = @(t) (1/(2*m)) * sum((X*t - y).^2);
gradFunction = @(t) (1/m) * X' * (X*t - y);

% fmincg 함수를 사용하여 theta(모델 파라미터) 값을 학습합니다.
[theta, cost] = fmincg(costFunction, initial_theta, options, gradFunction);


disp(theta);
