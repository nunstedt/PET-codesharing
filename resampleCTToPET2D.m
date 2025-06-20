function CTResampled = resampleCTToPET2D(CTIm, PETSize, CTFOV, PETFOV)
% RESAMPLECTTOPET - Crops and resamples a CT image to match PET resolution.
%
%   INPUT:
%       CTIm   - 3D matrix of the CT image (size: [CT_X, CT_Y])
%       PETSize - 2-element vector [PET_X, PET_Y] specifying PET image size (e.g., [256, 256])
%       CTFOV  - Field of view (FOV) of the CT image in mm (e.g., 500)
%       PETFOV - Field of view (FOV) of the PET image in mm (e.g., 300)
%
%   OUTPUT:
%       CTResampled - CT image cropped and resampled to match PET resolution
%

    % Compute voxel sizes
    CTVoxelSize = CTFOV / size(CTIm, 1);  % mm/pixel

    % Compute center indices for cropping
    centerIdx = floor(size(CTIm, 1) / 2);
    cropHalf = floor((PETFOV / 2) / CTVoxelSize);

    % Crop CT image around center
    CTCropped = CTIm(centerIdx - cropHalf : centerIdx + cropHalf - 1, ...
                         centerIdx - cropHalf : centerIdx + cropHalf - 1);
    
    % figure
    % imshow(CTCropped(:,:,25), [])

    % Resample cropped CT to match PET resolution
    CTResampled = imresize(CTCropped, [PETSize(1), PETSize(2)], 'bilinear');

end
