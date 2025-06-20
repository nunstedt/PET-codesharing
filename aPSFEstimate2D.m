function [kernelShape] = aPSFEstimate2D(coeffs,X,Y)
% APSFESTIMATE3D - Computes the 3D anisotropic point spread function (aPSF)
%                  using a sum of four exponentials.
%
%   INPUT:
%       coeffs - A vector of 8 parameters defining the exponential model:
%                [a, b, c, d, e, f, g, h], where each pair represents
%                the weight and decay factor of an exponential term.
%       X      - X-coordinates of the 2D grid (in mm).
%       Y      - Y-coordinates of the 2D grid (in mm).
%
%   OUTPUT:
%       kernelShape - Computed 2D kernel based on the given coefficients.


% Compute the 2D kernel using a sum of three exponentials.
kernelShape = coeffs(1) .* exp( sqrt(X.^2 + Y.^2) .* coeffs(2) ) + ...
              coeffs(3) .* exp( sqrt(X.^2 + Y.^2) .* coeffs(4) ) + ...
              coeffs(5) .* exp( sqrt(X.^2 + Y.^2) .* coeffs(6) );

% Compute the 2D kernel using a sum of four exponentials.
% kernelShape = coeffs(1) .* exp( sqrt(X.^2 + Y.^2) .* coeffs(2) ) + ...
%               coeffs(3) .* exp( sqrt(X.^2 + Y.^2) .* coeffs(4) ) + ...
%               coeffs(5) .* exp( sqrt(X.^2 + Y.^2) .* coeffs(6) ) + ...
%               coeffs(7) .* exp( sqrt(X.^2 + Y.^2) .* coeffs(8) );
end