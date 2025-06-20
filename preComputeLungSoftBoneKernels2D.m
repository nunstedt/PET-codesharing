function [lungKernel, softKernel, boneKernel] = preComputeLungSoftBoneKernels2D(pixelSize, kernelSize)
% PRECOMPUTELUNGSOFTBONEKERNELS - Computes positron range kernels for lung,
% soft tissue, and cortical bone using the sum of 4 exponentials as model 
% fitted to data.
%
%   INPUT:
%       pixelSize  - The physical size of each pixel (in mm).
%       kernelSize - The size of the kernel (number of pixels).
%
%   OUTPUT:
%       lungKernel - Precomputed kernel for lung tissue.
%       softKernel - Precomputed kernel for soft tissue.
%       boneKernel - Precomputed kernel for cortical bone.
%


% % Where model parameters come from
% F18
% load("dataF18/aPSF_models_Lung_Soft_CortF18.mat")
% lungCoeffs4Exp = coeffvalues(fLung_4Exp);
% softCoeffs4Exp = coeffvalues(fSoft_4Exp);
% boneCoeffs4Exp = coeffvalues(fCortical_4Exp);

% Exponential model parameters fitted from F18!! data in Carter et al
% lungCoeffs4Exp = 1.0e+02 .* [0.047512229545341  -1.189535152685521   0.004562271488979  -0.347836571142300   0.000324906595474  -0.091271943596803   0.000014071591954  -0.018629734712456];
% softCoeffs4Exp = 1.0e+02 .* [0.046599278930263  -3.776379516481885   0.004043852789245  -1.103362104617367   0.000319259572757  -0.297134075307484   0.000014299551444  -0.061687419372281];
% cortCoeffs4Exp = 1.0e+02 .* [0.054619619815922  -7.169918864823855   0.005270084803171  -2.205791965141861   0.000447137794947  -0.610672948461447   0.000017542033381  -0.115290989265292];
% 

% % Where model parameters come from
% Ga68
% load("dataGa68/aPSF_models_Lung_Soft_CortGa68.mat")
% lungCoeffs4Exp = coeffvalues(fLung_4Exp);
% softCoeffs4Exp = coeffvalues(fSoft_4Exp);
% boneCoeffs4Exp = coeffvalues(fCortical_4Exp);

% Exponential model parameters fitted from Ga68!! data in Carter et al
% lungCoeffs4Exp = [1.99289770145100	-14.0221025239714 -0.299167500510233 -2.16603309935822	0.443397451421166	-2.37896209787877	0.0143993393572844	-0.588901004671122];
% softCoeffs4Exp = [2.07653710515969	-45.3980265849574	0.00712423195570211	-8.74260369243776	0.101521415052282	-6.62658392548621	0.00793422756884768	-2.21221235994122];
% boneCoeffs4Exp = [2.02662334151737	-82.2175127126637	0.927916694176007	-14.8740183743410	-0.794158241460405	-13.9150474252753	0.0440056794404026	-5.88119292280254];

% 3exp was better R^2 fit!!!!!!!!!!!!!!!!!!!!!
lungCoeffs3Exp = [2.2817  -17.9492    0.3005   -4.1697    0.0154   -0.6053];
softCoeffs3Exp = [2.2764  -58.5195    0.3022  -13.6293    0.0154   -1.9790];
boneCoeffs3Exp = [2.4193 -108.3047    0.3225  -25.0723    0.0176   -3.8045];

% Compute the half size of the kernel for proper centering
half_size = floor(kernelSize / 2);

% Generate a coordinate grid centered at (0,0)
[X, Y] = meshgrid(-half_size:half_size, -half_size:half_size); 

% Convert indices to physical distances in mm
X = X .* pixelSize;
Y = flipud(Y) .* pixelSize; % Flipping Y to match MATLAB's coordinate convention (positive Y downward) Redundant due to circularly symmetric atm.

% Compute kernels using the fitted exponential model parameters.
lungKernelShape = aPSFEstimate2D(lungCoeffs3Exp,X,Y);
softKernelShape = aPSFEstimate2D(softCoeffs3Exp,X,Y);
boneKernelShape = aPSFEstimate2D(boneCoeffs3Exp,X,Y);

% Normalize the kernels so that their sum equals 1.
lungKernel = lungKernelShape ./ sum(lungKernelShape, 'all');
softKernel = softKernelShape ./ sum(softKernelShape, 'all');
boneKernel = boneKernelShape ./ sum(boneKernelShape, 'all');


end