function [kernelSize] = getPosRangeKernelSizeGa68(pixelSize,PSFThreshold)
% GETPOSRANGEKERNELSIZEGA68 - Computes the appropriate kernel size for the positron range
% blurring kernel in Ga-68 PET imaging, based on a given PSF threshold and voxel size.
%
%   INPUT:
%       pixelSize     - Size of a voxel edge in mm. Assumes cubic voxels (1D scalar).
%       PSFThreshold  - Threshold value at which the kernel should truncate the aPSF (e.g., 1e-6).
%
%   OUTPUT:
%       kernelSize    - The full kernel size in number of pixels (odd integer),
%                       chosen so that the PSF falls below the threshold.
%
%   NOTE:
%       Assumes 3D kernel (but actual dimensionality is not handled here),
%       using Carter et al. aPSF model for Ga-68 in lung tissue.


% Constants
max_r_Ga68 = 30; % [mm] - maximum radial extent to evaluate aPSF
r = linspace(10^-4, max_r_Ga68, 600); % [mm] - radial sampling (20 samples/mm)

% Coefficients for a 4-exponential fit of aPSF in lung tissue (from Carter et al.)
lungCoeffs4Exp = [1.99289770145100	-14.0221025239714 -0.299167500510233 -2.16603309935822 ...
                  0.443397451421166	-2.37896209787877	0.0143993393572844	-0.588901004671122];

% Function handle for aPSF estimate: sum of four exponentials
aPSFEst_r = @(r,coeffs) coeffs(1) .* exp(r .* coeffs(2)) + ...
                        coeffs(3) .* exp(r .* coeffs(4)) + ...
                        coeffs(5) .* exp(r .* coeffs(6)) + ...
                        coeffs(7) .* exp(r .* coeffs(8));

% Evaluate aPSF over the radial range
aPSFEstLung = aPSFEst_r(r, lungCoeffs4Exp);

% Find first index where aPSF drops below threshold
idx = find(aPSFEstLung < PSFThreshold, 1, 'first');
rThreshold = r(idx); % [mm] - radial distance where aPSF drops below threshold

% Compute number of pixels needed to cover up to rThreshold
half_nPixels = rThreshold / pixelSize;

% Round up to nearest even integer (ensures center pixel and symmetry)
half_nPixelsEven = 2 * ceil(half_nPixels / 2);

% Full kernel size is 1 (center) + 2 * half-width
kernelSize = 1 + 2 * half_nPixelsEven;

end