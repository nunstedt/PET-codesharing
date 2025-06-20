function [rhoIm] = CTIm2RhoEstIm(CTIm)
% CTIM2RHOESTIM Converts a CT image (Hounsfield Units) to an estimated 
%               density image (g/cm^3).
%   The conversion assumes the following bilinear relationship between CT 
%   numbers (HU) and density:
%     - CT = 1000 HU corresponds to a density of ~2 g/cm^3
%     - CT = 0 HU corresponds to a density of ~1 g/cm^3
%     - CT = -1000 HU corresponds to a density of ~1.293 × 10⁻³ g/cm^3 (air)
%
%
%   The transformation is modeled as two linear equations:
%     - For CT >= 0:    ρ = (1/1000) * CT + 1
%     - For CT < 0:     ρ = k1 * CT + m1, where k1 is derived from the (-1000, 1.293×10⁻³) to (0,1) mapping.


% Clips values in CT image to fit the domain the functions above are defined for.
CTIm(CTIm < -1000) = -1000;
CTIm(CTIm > 1000) = 1000;

% Compute slope and intercept for CT < 0
k1 = 9.98707e-04;  % (1 - 1.293e-3) / (0 - (-1000))
m1 = 1;            % Intercept at CT = 0

% Compute slope and intercept for CT >= 0
k2 = 10^-3; %1 / 1000;     % Slope from (0,1) to (1000,2)
m2 = 1;            % Intercept at CT = 0


% Define piecewise linear equations
y1 = k1 .* CTIm + m1;  % For CT < 0
y2 = k2 .* CTIm + m2;  % For CT >= 0


% Apply piecewise condition
rhoIm = y2;                 % Default: apply equation for CT >= 0
rhoIm(CTIm < 0) = y1(CTIm < 0);  % Apply equation for CT < 0

end