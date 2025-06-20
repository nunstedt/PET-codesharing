function [convdIm] = imPosRangeFilt_tripleImage2Dfaster(PETIm,padSize,lungKernel,softKernel,boneKernel,maskLungSoft,maskSoftBone,mixFactorLungSoft,mixFactorSoftBone)
% IMPOSRANGEFILT_TRIPLEIMAGE2Dfaster - Applies local positron range blurring to a 2D PET image
%                                      based on tissue density estimated from CT using the
%                                      triple image method.
%
%   INPUT:
%       PETIm         - 2D PET image matrix (activity image to be convolved).
%       padSize       - Size of padding applied to the image for boundary handling.
%       lungKernel    - Kernel for lung tissue convolution.
%       softKernel    - Kernel for soft tissue convolution.
%       boneKernel    - Kernel for bone tissue convolution.
%       maskLungSoft  - Mask for lung and soft tissue transition (logical array).
%       maskSoftBone  - Mask for soft and bone tissue transition (logical array).
%       mixFactorLungSoft - Mixing factor between lung and soft tissue (value between 0 and 1).
%       mixFactorSoftBone - Mixing factor between soft and bone tissue (value between 0 and 1).
%
%   OUTPUT:
%       convdIm       - Filtered PET image with local positron range blurring applied
%                       according to CT-derived tissue density.

% Replicate pad PETIm depending on size of kernel
% Apply padding using 'replicate' mode to extend edge values
paddedPETIm = padarray(PETIm, [padSize padSize], 'replicate', 'both');

% Computing 3 images (Triple image method)
convdImLung = convn(paddedPETIm, lungKernel, 'valid');
convdImSoft = convn(paddedPETIm, softKernel, 'valid');
convdImBone = convn(paddedPETIm, boneKernel, 'valid');


% Weighted combination of images - will this preserve sum of image?
% probably not? A convolution with a kernel summing to 1 would do that.
convdIm = maskLungSoft .* ((1 - mixFactorLungSoft) .* convdImLung + ... % if rho < 1 and mixfact (rho is close to 0) Lung dominates
                                mixFactorLungSoft .* convdImSoft)  + ... % and if rho is close to 1 Soft dominates
          maskSoftBone .* ((1 - mixFactorSoftBone) .* convdImSoft + ... % if rho > 1 and mask is close to 0 (rho is close to 1) Soft dominates 
                                mixFactorSoftBone .* convdImBone);       % and if mask is close to 1 (rho is close to 2) Bone dominates

end