% Fluid properties:
V = 10;  % Fluid speed
rho = 1.225;  % Fluid density
mu = 16.82e-6;  % Dynamic viscosity

% Desired/assumed turbine properties:
P = 1e3; % Power output
Cp = 0.3 * 0.75; % Overall efficiency (blades*generator)
tsr = 5;  % Tip speed ratio  %TODO determine based on sound
B = 3;  % Number of blades

% Hub geometry:
hub_radius = 0.20;  % Radius at the end of the shaft
shaft_length = 0.1;
shaft_diameter = 0.06;
trans_length = 0.1;  % Length of transition between shaft and blade
fat_length = 0.3;  % Length of transition from fat profile to main profile
tip_length = 0;  % Length of tip (section where better airfoil will be introduced)
resolution = 120;  % Number of segments (lower for faster simulation, higer for final print)

% Create a databse of possible foils (see 'f' function for property order):
foils = struct(...
    'circle', f(5, 0, 0.47, 1, 0), ...
    'NACA0018', f(5, 0.5215, 0.01023, 0.18, 3), ...
    'NACA4418', f(2.5, 0.7263, 0.00945, 0.18, 2), ...
    'NACA4418_best_angle', f(6.5, 1.1326, 0.01214, 0.18, 2));

% Foils to actually use (in order):
FS = [foils.circle, foils.NACA0018, foils.NACA4418, foils.NACA4418_best_angle];




% Radius set based on desired power:
R = sqrt((2*P)/(Cp*rho*pi*V^3))  % Total blade-tip radius

% Calculated turbine properties:
w = (tsr*V/R);  % Angular velocity (at desired tsr)
rpm = 60 * w/(2*pi)  % Convert angular velocity to rpm
T = P/w  % Required torque
%fprintf('Tip speed: %f km/h\n\n', (tsr*V)*60*60/1000);  % Indicative of turbine loudness
fprintf('Tip speed: %f m/s\n\n', tsr*V);


% Split up blade into sections (for BEM theory):
rs = linspace(hub_radius-shaft_length, R, resolution);  % Radii
rs(end) = rs(end)*0.9995;  % Avoid NaN errors at tip

% Assign an airfoil index to each blade element (decimal index allowed, causing interpolation):
[~, trans_i] = min(abs(rs - (hub_radius)));
[~, fat_i] = min(abs(rs - (hub_radius+trans_length)));
[~, main_i] = min(abs(rs - (hub_radius+trans_length+fat_length)));
[~, tip_i] = min(abs(rs - (R-tip_length)));
fs = [ones(1,trans_i), sinspace(1,2,fat_i-trans_i), sinspace(2,3, main_i-fat_i), ones(1,tip_i-main_i)*3, sinspace(3,4, resolution-tip_i)];
%fs = [ones(1,trans_i), linspace(1,2,fat_i-trans_i), linspace(2,3, main_i-fat_i), ones(1,resolution-main_i)*3];  % Without tip section

rs(trans_i) = hub_radius;  % Ensures shaft length is exact irrespective of resolution


% Calculate interpolated foil properties:
bs = prop_lerp(fs, FS, 'b');
Cls = prop_lerp(fs, FS, 'Cl');
Cds = prop_lerp(fs, FS, 'Cd');
ths = prop_lerp(fs, FS, 'th');



 %%%  Determine optimal blade geometry (iteratively):  %%%

% Initially, induction factors are assumed to be zero. These are adjusted
% iteratively.
ind_as = zeros(size(rs));
ind_rs = ind_as;

% Main iteration loop:
for i = 1:1000
    ts = twist(rs, V, w, bs, ind_as, ind_rs);  % Calculate twist from relative wind velocity and angle of attack
    cs = chord_length(rs, ts-bs, B, tsr, R);  % Optimal chord-length is approximated by formula
    
    sols = (B*cs)./(2*pi*rs);  % Local 'blade solidities' (Dimensionless, used for calculations)
    
    % Save current induction factors for later comparison (to check convergence):
    old_a = ind_as;
    old_r = ind_rs;
    
    % Improve induction factor estimates given new twists and chord lengths:
    Q = tip_loss_correction(rs, B, ts-bs, R);  % Tip loss coefficient (rounds off blade tip)
    ind_as = a_induction_factor(sols, ts-bs, Q, Cls, Cds);  % Axial induction factor
    ind_rs = r_induction_factor(sols, ts-bs, Q, Cls, Cds, ind_as, rs, tsr, R);  % Radial induction factor
    %mean_a = mean(ind_a)  % For debugging
    %mean_r = mean(ind_r)
    
    % Check whether we've converged:
    if sum(abs(ind_as-old_a))+sum(abs(ind_rs-old_r)) <= 1e-6
        fprintf('Solution converged in %i iterations.\n', i);
        break;
    elseif i==1000
        fprintf('Exceeded maximum iterations. Did not converge!\n');
    end
end

% Generate shaft geometry:
cs(1:trans_i) = shaft_diameter;
cs(trans_i:fat_i) = lerp(shaft_diameter, cs(fat_i), linspace(0,1,fat_i-trans_i+1) );

% Calculate Reynolds numbers (useful for manual profile selection):
Reys = rho*local_tsr(rs, tsr, R)*V.*cs ./ mu;

%V = 10  % Test blade performance at different wind speed

% Calculate torque and drag using numerical integration:
dT = @(n) delta_torque(rs(n), ts(n)-bs, sols(n), ind_as(n), V, rho, Cls, Cds);
dTs = dT(1:length(rs));
torque = trapez_int(rs, dTs);

dD = @(n) delta_drag(rs(n), ts(n)-bs, sols(n), ind_as(n), V, rho, Cls, Cds);
dDs = dD(1:length(rs));
drag = trapez_int(rs, dDs);

% Very rough stress calculation:
c_ts = cs .* ths * 0.5;  % Mean chord thickness (estimate!)
th = 3e-3;  % Shell thickness
dens = 1.11e3;  % Material density
c_As = (2*cs + 2*c_ts)*th;  % Cross sectional areas

M_ds = rs * 0;  % Moments due to drag
M_ts = M_ds;  % Moments due to torque
S_cs = M_ds;  % Tension stresses due to centripital force  %TODO
for i = 1:length(rs)
    M_ds(i) = -trapez_int(rs(end:-1:i), rs(end:-1:i) .* dDs(end:-1:i)/B);
    M_ts(i) = -trapez_int(rs(end:-1:i), dTs(end:-1:i)/B);
end
cent_stresses = w^2 * (rs(1:end-1)+rs(2:end))/2 .* (rs(2:end)-rs(1:end-1)) * dens;  % Centripital force on each blade element
for i = 1:length(rs)-1
    S_cs(i) =  sum(cent_stresses(end:-1:i));
end

s = stress(cs, c_ts, M_ds, M_ts, S_cs, th);

% Done! Print results:

% Plot performance properties:
hold off;
subplot(4,1,1);
plot(rs, s*1e-6);
ylabel('Stress (MPa)');
subplot(4,1,2);
plot(rs, dTs);
ylabel('Lift (N/m)');
subplot(4,1,3);
plot(rs, dDs);
ylabel('Drag (N/m)');
subplot(4,1,4);
plot(rs, Reys);
ylabel('Reynold''s number');
xlabel('Position along blade (m)');

% Export csv files:
%  Blender:
blade_data = [rs-rs(1); cs; 90-rad2deg(ts); fs].';
csvwrite('blade_blend.csv', blade_data);
%  QBlade (does not support interpolation):
blade_data = [rs-rs(1); cs; 90-rad2deg(ts); FS(round(fs)).('i')].';
csvwrite('blade_qblade.csv', blade_data);

torque
drag

fprintf('Hub radius: %fm\n', hub_radius);
fprintf('Inner radius: %fm\n', hub_radius-shaft_length);

fprintf('\nRe:\n Mean:%0.0f, Min:%0.0f, Max:%0.0f\n\n', mean(Reys), min(Reys), max(Reys));

fprintf('Maximum stress (estimate): %f MPa\n\n', max(abs(s)) * 1e-6);

final_efficiency = (2*torque*w)/(pi*R^2*rho*V^3)
final_power = (pi*R^2*rho*V^3 * final_efficiency)/2

 %%%  Functions  %%%

% Speed ratio at any radius:
function l = local_tsr(r, tsr, R)
    l = r*tsr/R;
end

% Blade twist at radius:
function a = twist(r, V, w, b, ind_a, ind_r)
    a = atan( (w*r.*(1+ind_r)) ./ (V*(1-ind_a)) ) + b;  % From BME paper
end

% Chord length at radius:
function c = chord_length(r, phi, B, tsr, R)
    l = local_tsr(r,tsr,R);
    c = (8*pi.*r.*cos(phi))./(3*B.*l);  % From BME paper
end

% Tip loss correction factor:
function Q = tip_loss_correction(r, B, phi, R)
    frac = -B/2 .* (1-r/R) ./ ((r/R).*cos(phi));
    Q = 2/pi .* acos(exp(frac));
end

% Axial induction factor:
function ind_a = a_induction_factor(s, phi, Q, Cl, Cd)  %s is solidity, Q is tip loss correction factor
    frac = s.*(Cl.*sin(phi)+Cd.*cos(phi)) ./ (4*Q.*cos(phi).^2);
    ind_a = frac ./ (1+frac);
end

% Radial induction factor:
function ind_r = r_induction_factor(s, phi, Q, Cl, Cd, ind_a, r, tsr, R)
    ltsr = local_tsr(r, tsr, R);
    %frac = s.*(Cl.*sin(phi)-Cd.*cos(phi)) ./ (4*Q.*ltsr.*cos(phi).^2);  %old typo
    frac = s.*(Cl.*cos(phi)-Cd.*sin(phi)) ./ (4*Q.*ltsr.*cos(phi).^2);
    ind_r = frac .* (1-ind_a);
end

% Torque on blade element
function dT = delta_torque(r, phi, sol, ind_a, V, rho, Cl, Cd)
    dT = sol .* pi*rho .* (V^2*(1-ind_a).^2)./cos(phi).^2 .* (Cl.*cos(phi)-Cd.*sin(phi)).*r.^2;
end

% Drag on blade element:
function dD = delta_drag(r, phi, sol, ind_a, V, rho, Cl, Cd)
    dD = sol .* pi*rho .* (V^2*(1-ind_a).^2)./cos(phi).^2 .* (Cl.*sin(phi)+Cd.*cos(phi)).*r;
end

% Peak stress at blade element:
function s = stress(c, t, M_d, M_t, S_cs, th)
    I_d = (t.^2 .* c * th) / 3;
    I_t = (c.^2 .* t * th) / 3;
    s = (M_d .* t/2)./I_d + (M_t .* c/2)./I_t + S_cs;
end

% Creates a struct out of foil data:
function out = f(b, Cl, Cd, th, i)
    % Order of blade properties:
    %  Angle of attack, Lift coefficient, Drag coefficient, max thickness, qBlade foil index
    out = struct('b',deg2rad(b), 'Cl',Cl, 'Cd',Cd, 'th',th, 'i',i);
end

% Linearly interpolate values:
function x = lerp(a, b, f)
    x = a.*(1-f) + b.*f;
end

% Calculate properties of an interpolated foil:
function x = prop_lerp(f, F, name)
    x = lerp( [F(floor(f)).(name)], [F(ceil(f)).(name)], mod(f,1) );  % Linear
end

% Like 'linspace', but sinusoidal instead of linear (for smoother transitions):
function x = sinspace(a, b, n)
    x = a+(1-(cos( linspace(0,1,n) *pi)+1)/2)*(b-a);
end
