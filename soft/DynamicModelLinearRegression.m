
% git passwd = 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
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

% i think it equals to
% Y = aX_{1} + bX_{1}X_{2} + cX_{1}*X_{1}