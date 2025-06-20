function [padSize, lungKernel, softKernel, boneKernel, maskLungSoft, maskSoftBone, mixFactorLungSoft, mixFactorSoftBone] = preComputeAllPosRangeKernelParams(PSFThreshold,pixelSize,CT2D_25_resampled)
    
    % Estimate tissue density from CT values
    rhoIm = CTIm2RhoEstIm(CT2D_25_resampled);
    % figure, imshow(rhoIm,[])
    % Mask for lung and soft tissue transition
    maskLungSoft = (rhoIm < 1); % Logical array: 1 where rho < 1, 0 otherwise
    % Mask for soft and bone tissue transition
    maskSoftBone = (rhoIm >= 1); % Logical array: 1 where rho ≥ 1, 0 otherwise
    
    % Tissue with density between lung and soft tissue (rho in [0,1))
    mixFactorLungSoft = rhoIm;
    mixFactorLungSoft(mixFactorLungSoft >= 1) = 0; % zero for all elements where rho > 1
    % Tissue with density between soft and bone tissue (rho in [1,2])
    mixFactorSoftBone = rhoIm - 1; % values between 0 and 1
    mixFactorSoftBone(mixFactorSoftBone < 0) = 0; % zero for all elements where rho < 0
    
    % Determine kernel size so that the positron range kernel includes all values where aPSF ≥ PSFThreshold
    kernelSize = getPosRangeKernelSizeGa68(pixelSize, PSFThreshold); % odd
    fprintf('\nPositron kernel size: %d\n', kernelSize)
    
    % Compute padding size
    padSize = floor(kernelSize / 2);
    
    % Precompute 2D positron range kernels for lung, soft tissue, and bone. 
    [lungKernel, softKernel, boneKernel] = preComputeLungSoftBoneKernels2D(pixelSize, kernelSize);



end