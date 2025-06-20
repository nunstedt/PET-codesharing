clear all, close all, imtool close all
clc


% Adding paths
% addpath('simulation')

%%
sigmaConst = (1/sqrt(8*log(2)));          % 1/2.3548
cSpeed     = physconst('LightSpeed')*1000;      % Speed of light mm/sec

%% 1) Comment this!
gNum       = 1E3;
sysRad     = 125;                % mm
sysCTR     = 100E-12;            % ps
fwhmDetPE  = 1;                  % mm
sysFOV     = 300;
boxFactor  = 2.5;                % xFWHM
matRad     = 4*32;
matSize    = 2*matRad+1;
mm2vox     = matSize/(sysFOV);
sensAngles = floor(matSize/2); %from 3d sim floor(pi*sysRad/max(fwhmDetPE,4*sysRad/matSize));
load('PETCT_2D_Huber_25.mat');


%% And this!
fwhmTmm     = cSpeed*sysCTR/2;    % mm
sigmaTmm    = fwhmTmm*sigmaConst; % mm
sigmaDetmm  = fwhmDetPE*sigmaConst;      % mm
fwhmAcol    = 0.5*pi/180;         % rad
sysFWHMmm   = sqrt((sysRad*tan(fwhmAcol/2))^2 + (1.5*sigmaDetmm/2/sigmaConst)^2); 
sigmaXmm    = sysFWHMmm*sigmaConst;
varMatmm    = [ sigmaXmm^2   0 ; 0 sigmaTmm^2   ];
iVarMatmm   = [ 1/sigmaXmm^2 0 ; 0 1/sigmaTmm^2 ];
varDetmm    = det(varMatmm);
normSplatmm = 1/(2*pi*sqrt(varDetmm));

boxSizemm  = boxFactor*max(sysFWHMmm,fwhmTmm)*[1 1];
boxSize    = ceil(matSize/sysRad * boxSizemm);
boxSize    = 2*floor(boxSize/2) + mod(boxSize,2);
boxStride  = int32(boxSize(1));





%% Grid 
xyGridmm    = linspace(-sysFOV/2,sysFOV/2,matSize);
yxGridmm    = linspace(-sysFOV/2,sysFOV/2,matSize);
[YXmm,XYmm] = meshgrid(xyGridmm,xyGridmm); % note x/y swap
YXmm        = single(YXmm);
XYmm        = single(XYmm);
supDisk   = single(xyGridmm(ones(matSize(1),1),:).^2 + xyGridmm(ones(1,matSize(1)),:)'.^2 <= floor(sysRad)^2);
sup1      = single(xyGridmm(ones(matSize(1),1),:).^2 + xyGridmm(ones(1,matSize(1)),:)'.^2 <= floor(sysRad-5)^2);
nSup1 = sum(sup1(:));


%% Centering PET and CT images
% Correctly shifted matSize x matSize ground truth
phantPETGT = imtranslate(PET2D_Huber_25, [9,-7]);
phantPETGT = imresize(phantPETGT,[matSize, matSize],"bilinear");


minvalabs = abs(min(CT2D_25(:)));
CT2D_25_tmp = CT2D_25 + minvalabs;
CT2D_25_shifted = imtranslate(CT2D_25_tmp, [10.8, -8.4]);
CT2D_25_shifted = CT2D_25_shifted - minvalabs;


%% Build simulation phantom

phantBoost = 3;
boostedMatSize = phantBoost*matSize;

phantPETSim  = imresize(phantPETGT,[boostedMatSize, boostedMatSize],'bilinear');

boostedPhantSize = size(phantPETSim);
xyGridBoost    = linspace(-sysFOV/2,sysFOV/2,boostedPhantSize(1));
yxGridBoost    = linspace(-sysFOV/2,sysFOV/2,boostedPhantSize(2));
supBoost       = single(xyGridBoost(ones(boostedPhantSize(1),1),:).^2 + xyGridBoost(ones(1,boostedPhantSize(1)),:)'.^2 <= ...
    floor(sysRad-5)^2);

% HAVE THIS HERE TOO?
% phantPETSim = supBoost.*phantPETSim;
% phantPETSim(phantPETSim<max(phantPETSim(:))/500) = 0;



%% Cropping and resampling CT image for positron range and attenuation for Simulation!!!

fprintf('\nCropping and resampling CT image for attenuation in simulation...\n')
tic
PETFOV = 300; % matSize/mm2vox; % [mm].
CTFOV = 500; % [mm]
phantPETSizeSim = size(phantPETSim);
pixelSizeSim = PETFOV / boostedPhantSize(1);
CTenergy = 70; % ?????????


% Crop and resample CT image to PET dim. and size.
CT2D_25_resampledSim = resampleCTToPET2D(CT2D_25_shifted,[boostedPhantSize(1), boostedPhantSize(2)], CTFOV, PETFOV); 

% Create mu map for simulation
muCTSim = single(ct2mu(CT2D_25_resampledSim, CTenergy)); 
toc

%% Positron range

fprintf('\nApplying positron range kernel to upsampled ground truth phantom for simulation...')

% Threshold value at which the kernel should truncate aPSF (e.g., 1e-6).
PSFThresholdSim = 10^-4; 
tic
% Precompute inputs for imPosRangeFilt_tripleImage2Dfaster for simulation
[padSizeSim, lungKernelSim, softKernelSim,...
 boneKernelSim, maskLungSoftSim, maskSoftBoneSim,...
 mixFactorLungSoftSim, mixFactorSoftBoneSim] = ...
 preComputeAllPosRangeKernelParams(PSFThresholdSim,pixelSizeSim,CT2D_25_resampledSim);

% % FOR PLOTS!
% phantPETSim_norm = phantPETSim / max(phantPETSim(:));

% Apply positron range conv. to ground truth before simulation
phantPETSim_PR = imPosRangeFilt_tripleImage2Dfaster(phantPETSim, padSizeSim,...
            lungKernelSim, softKernelSim, boneKernelSim, maskLungSoftSim, maskSoftBoneSim,...
            mixFactorLungSoftSim, mixFactorSoftBoneSim);

% ?!?!?!!!!!DO BOTH OF THESE HERE OR NO????!!!!!!!!!
phantPETSim_PR = supBoost.*phantPETSim_PR;
phantPETSim_PR(phantPETSim_PR<max(phantPETSim_PR(:))/500) = 0; % removes
% some effect from PR?!?
toc

%%

% Create output folder for saving plots
% outputFolder = fullfile(pwd, 'report_results','resultplotsRich','positron_range_sim');
% if ~exist(outputFolder, 'dir')
%     mkdir(outputFolder);
% end
% 
% % ----------- Normalized GT PET image -----------
% fig = figure('Color', 'white');
% imshow(phantPETSim_norm, []);
% colormap(flipud(gray));
% colorbar;
% axis off tight;  % Removes ticks but keeps room for colorbar
% 
% filePath = fullfile(outputFolder, 'upsampled_GT_PET_norm.png');
% print(fig, filePath, '-dpng', '-r300');
% close(fig);
% 
% % ----------- Positron-ranged version -----------
% fig = figure('Color', 'white');
% imshow(phantPETSim_PR, []);
% colormap(flipud(gray));
% colorbar;
% axis off tight;  % Removes ticks but keeps room for colorbar
% 
% filePath = fullfile(outputFolder, 'pos_ranged_upsampled_GT_PET_norm.png');
% print(fig, filePath, '-dpng', '-r300');
% close(fig);
% 
% % ----------- Difference image -----------
% diff_im = phantPETSim_norm - phantPETSim_PR;
% 
% fig = figure('Color', 'white');
% imshow(diff_im, []);
% colormap(flipud(gray));
% colorbar;
% axis off tight;  % Removes ticks but keeps room for colorbar
% 
% 
% filePath = fullfile(outputFolder, 'diffIm_upsampld_GT_minus_upsampled_GT_PET_pos_ranged.png');
% print(fig, filePath, '-dpng', '-r300');
% close(fig);

%% For positron range simulation results

% % OLD!!!!!!!!!!!!!!!
% 
% % Create output folder for saving plots
% outputFolder = fullfile(pwd, 'report_results','resultplotsRich','positron_range_sim');
% if ~exist(outputFolder, 'dir')
%     mkdir(outputFolder);
% end
% 
% figure,imshow(phantPETSim_norm,[])
% filePath = fullfile(outputFolder, 'upsampled_GT_PET_norm.png');
% colormap(flip(gray));
% colorbar;
% exportgraphics(gcf, filePath);
% close(gcf);
% 
% 
% 
% figure,imshow(phantPETSim_PR,[])
% filePath = fullfile(outputFolder, 'pos_ranged_upsampled_GT_PET_norm.png');
% colormap(flip(gray));
% colorbar;
% exportgraphics(gcf, filePath);
% close(gcf);
% 
% 
% 
% % figure, imshow(phantPETSim,[])
% diff_im = (phantPETSim_norm-phantPETSim_PR);
% % figure,imshow(phantPETSim,[])
% figure, imshow((diff_im),[])
% colormap(flip(gray));
% colorbar;
% 
% filePath = fullfile(outputFolder, 'diffIm_upsampld_GT_minus_upsampled_GT_PET_pos_ranged.png');
% exportgraphics(gcf, filePath);
% close(gcf);


%% Simulation pixel-angle-method

sensSimAnglesFactor = 3; %ross said 3-10x on slack?   % Factor finer sampled angles for Simulation than sensitivity image
numSimAngles = sensSimAnglesFactor .* sensAngles; %  % Total number of angles sampled between -pi/2 and pi/2 - same as sensitivity image angles for now

% % % % % TEMPORARY
% % % % % % % % numSimAngles = 18;

% fprintf('\nRunning generateDataVecSim_Attn()... ')
% [X0mmSim, X1mmSim, X2mmSim, X0rlmmSim, ThetasSim, dataVecSim, muVecSim, idxsNonZeroGT] = ...
%     generateDataVecSim_Attn(phantPETSim_PR, xyGridBoost, yxGridBoost, ...
%     muCTSim, PETFOV, matSize, sysRad, numSimAngles);

fprintf('\nLoading saved generateDataVecSim_Attn() data... ')
tic
load("report_results\simData_pixelanglemethod\sim_pix_angl_method_3xim3xangls_r_125_sysfov300.mat")
toc
% save("report_results\simData_pixelanglemethod\sim_pix_angl_method_3xim3xangls_r_125_sysfov300.mat", "X0mmSim", "X1mmSim", "X2mmSim", "X0rlmmSim", "ThetasSim",...
%     "dataVecSim", "muVecSim", "idxsNonZeroGT", "phantPETSim_PR", "xyGridBoost", "yxGridBoost", ...
%     "muCTSim", "PETFOV", "matSize", "sysRad", "numSimAngles", '-v7.3')

%%

fprintf('\nRunning sampleDataVecSim_Attn()... ')
[dataVecSampled, X0mmSampled, X1mmSampled, X2mmSampled, X0rlmmSampled, ThetasSampled, muVecSampled, numG] = ...
    sampleDataVecSim_Attn(dataVecSim, X0mmSim, X1mmSim, X2mmSim, ...
    X0rlmmSim, ThetasSim, muVecSim, phantPETSim_PR, idxsNonZeroGT, gNum, sysRad, numSimAngles);

% save("report_results\simData_pixelanglemethod\sampl_pix_angl_method_10_7counts_from_3xim3xangls.mat", "dataVecSampled", "X0mmSampled", ...
    % "X1mmSampled", "X2mmSampled", "X0rlmmSampled", "ThetasSampled", "muVecSampled", "numG", "gNum", "sysRad", "numSimAngles", '-v7.3')

% load("report_results\simData_pixelanglemethod\sampl_pix_angl_method_10_7counts_from_3xim3xangls.mat")





%% Resolution loss mechanisms

fprintf('\nApplying resolution loss mechanisms... \n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD 2D RES LOSS MECHANISMS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% non-collinearity 2D

tic
[X1C,X2C,X0C] = nonColl2D(X1mmSampled',X2mmSampled',X0mmSampled',sysRad);
nonColl2DTime = toc;
fprintf("2D non-collinearity took    %011.6f seconds.\n", nonColl2DTime)

% % RICHARD PLOTTING STUFF
% if numG<1E5
%     figure,
%     plot(X1C(:,1),X1C(:,2),'.g'),hold on
%     plot(X2C(:,1),X2C(:,2),'.r'),axis square
% 
%     plot(X0C(:,1), X0C(:,2), '.b')     % as this one
% 
%     xlim(1.1*sysRad*[-1 1])
%     ylim(1.1*sysRad*[-1 1])
%     axis equal
%     title("X0,X1,X2 after non-coll")
% end



% detector resolution loss 2D
tic
[X1D,X2D,X0D,tC] = detResLoss2D(X1C,X2C,X0C,sysRad,sigmaDetmm);
detResLossTime = toc;
fprintf("2D detector res loss took   %011.6f seconds.\n", detResLossTime)


% % RICHARD PLOTTING STUFF
% if numG<1E5
%     figure,
%     plot(X1D(:,1),X1D(:,2),'.g'),hold on
%     plot(X2D(:,1),X2D(:,2),'.r'),axis square
% 
%     % plot(dataVec(2,:),dataVec(3,:),'.y') % should be same
%     plot(X0D(:,1), X0D(:,2), '.b')     % as this one
% 
%     xlim(1.1*sysRad*[-1 1])
%     ylim(1.1*sysRad*[-1 1])
%     axis equal
%     title("X0,X1,X2 after detResLoss")
% end



% CTR res loss 2D
tic
[dT,X0F]    = tofResLoss2D(X1D,X2D,tC,sigmaTmm);
tofResLossTime = toc;
fprintf("2D TOF res loss took        %011.6f seconds.\n", tofResLossTime)


% % RICHARD PLOTTING STUFF
% if numG<1E5
%     figure,
%     plot(X1D(:,1),X1D(:,2),'.g'),hold on
%     plot(X2D(:,1),X2D(:,2),'.r'),axis square
% 
%     % plot(dataVec(2,:),dataVec(3,:),'.y') % should be same
%     plot(X0F(:,1), X0F(:,2), '.b')     % as this one
% 
%     xlim(1.1*sysRad*[-1 1])
%     ylim(1.1*sysRad*[-1 1])
%     axis equal
%     title("X0,X1,X2 after resolution loss")
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some more data processing %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% REMOVE SIM POINTS OUTSIDE AND CLOSE TO RING!
indOutside = checkDetectorRadii(X1D,X2D,X0F,sysRad);
X0F(indOutside,:) = 0.5.*X1D(indOutside,:) + 0.5.*X2D(indOutside,:);
indOutside = checkDetectorRadii(X1D,X2D,X0F,sysRad);
X0F(indOutside,:) = 0;
indOutside = checkDetectorRadii(X1D,X2D,X0F,sysRad);


% flip detectors so that r always closest to X1
indXLR      = find( vecnorm(X0F-X1D,2,2) > 0.5*vecnorm(X2D-X1D,2,2));
X1tmp       = zeros(numG,2);
X2tmp       = zeros(numG,2);
X1tmp(indXLR,:) = X2D(indXLR,:);
X2tmp(indXLR,:) = X1D(indXLR,:);
X1D(indXLR,:) = X1tmp(indXLR,:);
X2D(indXLR,:) = X2tmp(indXLR,:);


% randomize order of data %
numG    = numel(tC);
indxT   = randperm(numG);
X1T     = X1D(indxT,:);
X2T     = X2D(indxT,:);
X0T     = X0F(indxT,:);
dT      = dT(indxT);
% deltaT = dT(indxT)*voxSizeS(1);
% tAC    = tAC(indxT); % attenuation coeff from sim, dont need?

% get r, l, theta %
r       = vecnorm(X0T-X1T,2,2);
l       = vecnorm(X2T-X1T,2,2);
thetas  = atan2(X2T(:,2)-X1T(:,2),X2T(:,1)-X1T(:,1));



fprintf('\nSimulation completed! \n')

% X0s FROM BRAIN IS ORIENTED CORRECTLY WHEN USING PLOT IN OUTPUT HERE!!!!


% --------------- PART OF SIMULATION END -------------------
%% Cropping and resampling CT image for positron range and attenuation for 
%  IMAGE RECONSTRUCTION!!!

fprintf('\nCropping and resampling CT image for attenuation in simulation...\n')
% % Variables needed for CT resampling and crop.
PETFOV = 300; % = matSize/mm2vox; % [mm].
CTFOV = 500; % [mm]
CTenergy = 70; % ?????????
pixelSize = PETFOV/matSize; % [mm/pixel]

% Crop and resample CT image
CT2D_25_resampled = resampleCTToPET2D(CT2D_25_shifted,[matSize, matSize], CTFOV, PETFOV); % 257x257 = PET size

% Create mu map
muCT = single(ct2mu(CT2D_25_resampled, CTenergy)); % args CTim, CTenergy

% figure, imshowpair(phantPETGT, CT2D_25_resampled, 'blend')

%% Positron range for image recon.

fprintf('\nPrecomputing positron range kernel parameters for image recon...')
% Threshold value at which the kernel should truncate the aPSF (e.g., 1e-6).
PSFThreshold = 10^-4; 

% Precompute inputs for imPosRangeFilt_tripleImage2Dfaster
[padSize, lungKernel, softKernel,...
boneKernel, maskLungSoft, maskSoftBone,...
mixFactorLungSoft, mixFactorSoftBone] = preComputeAllPosRangeKernelParams(PSFThreshold, pixelSize, CT2D_25_resampled);

%%

fprintf('\nCalculate attenuation for LORs in dataVec...\n')
tic
% Get attenuation values exp(-int mu(s) ds) for all events
muVec = getAttenuationMEX(X1T, X2T,single(PETFOV),single(PETFOV),single(3.75),flipud(muCT),int32(matSize),int32(matSize),int32(numG)); % ONLY FLIPUD ON muCT SINCE WE X1,X2 IS CORRECTLY ORIENTED, ONLY X0 WAS FLIPPED FOR RECONSTRUCTION!!!!
toc


%% Precomputing covariance matrices

% process data to fit brainSimfast_2D

% FLIPPING TO ORIENT TO ROSS'S WAY FOR RECON
X1mm        = flipud(X1T'); 
X1mm(1,:)   = -X1mm(1,:);
X2mm        = flipud(X2T'); 
X2mm(1,:)   = -X2mm(1,:);
X0mm        = flipud(X0T'); 
X0mm(1,:)   = -X0mm(1,:);


% % ROSS PLOTTING STUFF
% if numG<1E5
%     figure,
% % % % %     plot(X1mm(2,:),-X1mm(1,:),'.g'),hold on X1,X2 ARE NOT FLIPPED ABOVE SO WILL SHOW INCORRECTLY!
% % % % %     plot(X2mm(2,:),-X2mm(1,:),'.r'),axis square
%     xlim(1.*sysRad*[-1 1])
%     ylim(1.*sysRad*[-1 1])
%     plot(X0mm(2,:),-X0mm(1,:),'.b')
%     title("Ross code")
% 
% end

%%

fprintf('\nPrecompute covariance matrices... \n')
[dataVec, dataVecG, varMatLImm, iVarMatLImm, varNormLImm, meanVarLImm] = ...
    preComputeCovarianceMatrices(numG, X0mm, thetas-pi/2, r, l, ... %%% -pi/2!!!!!!!
    sigmaConst, fwhmAcol, sigmaDetmm, sigmaTmm, ...
    varMatmm, iVarMatmm, normSplatmm);

% %% LI Box AVX2 mex forward test
% fprintf('\tLocally  invariant forward-projector (Box+Chunk+AVX2)  ... ')
% imgIn = supDisk;
% tic
% forTmp = forwardProjLInv_BoxChunk(dataVec,imgIn, ...
%     single(xyGridmm),single(yxGridmm),single(mm2vox),boxStride,int32(32));
% forwardDataLi_xBoxC_c = dataVec;
% forwardDataLi_xBoxC_c(1,:) = forTmp;
% toc
% 
% %% LI Box AVX2 mex test
% fprintf('\tLocally  invariant back-projector    (Box+Chunk+AVX2)  ... ')
% dataIn = forwardDataLi_xBoxC_c;
% tic
% backImgLi_xBoxC_c = backProjLInv_BoxChunk(dataIn, ...
%     single(xyGridmm),single(yxGridmm),single(mm2vox),boxStride,int32(32));
% backImgLi_xBoxC_c = supDisk.*backImgLi_xBoxC_c;
% toc
% % %% Compute comparison acclerations and norms
% 
% IMG = double(backImgLi_xBoxC_c);
% find(IMG)
% imshow(IMG/max(IMG(:)),[])
% imtool(1-IMG/max(IMG(:)),[]);

% optVoxmm        = sysFWHMmm/3;
% fwhmXY          = sqrt(meanVarLImm(1,1))/sigmaConst;
% fwhmXYrecov     = max(fwhmXY,2/mm2vox);
% fwhmCTR         = sqrt(meanVarLImm(2,2))/sigmaConst;
% SNRvolml        = 4/3*pi*fwhmXY^2*fwhmCTR/8/1000;
% SNRvol20ml      = 4/3*pi*fwhmXY^2*min(fwhmCTR,200)/8/1000;
% 
% SNRvolmlrecov   = 4/3*pi*fwhmXYrecov^2*fwhmCTR/8/1000;
% SNRvol20mlrecov = 4/3*pi*fwhmXYrecov^2*min(fwhmCTR,200)/8/1000;
% 
% fprintf('\nSimulation parameters:\n')
% fprintf('\tDetection events = %7.1e\n',numG)
% fprintf('Scanner model:\n')
% fprintf('\tRadius           = %6.1f mm\n',sysRad)
% fprintf('\tCTR              = %6.1f ps\n',sysCTR*1E12)
% fprintf('\tFWHM(t)          = %6.2f mm\n',fwhmCTR)
% fprintf('\tFWHM(x)          = %6.2f mm\n',fwhmXY)
% fprintf('\tEllipsoid Vols.  = %7.4f (%7.4f) ml\n',SNRvolml,SNRvol20ml)
% 
% fprintf('Reconstruction model:\n');
% fprintf('\tBox size         = (%6.1f,%6.1f) mm\n',boxSizemm(1),boxSizemm(2))
% fprintf('\tVoxel size       = (%6.3f,%6.3f) mm\n',1/mm2vox,1/mm2vox)
% fprintf('\tVoxel opt        = (%6.3f,%6.3f) mm\n', optVoxmm,optVoxmm)
% fprintf('\tBox size         = (%4d,%4d)     vox\n',boxSize(1),boxSize(2))
% fprintf('\tMatrix size      = (%4d,%4d)     vox\n',matSize,matSize)
% fprintf('\tMatrix opt       = (%4d,%4d)     vox\n',round((2*sysRad+1)/optVoxmm),round((2*sysRad+1)/optVoxmm))
% fprintf('\tFWHM recov.      = %6.2f          mm\n',fwhmXYrecov)
% fprintf('\tEllipsoid recov. = %7.4f (%7.4f) ml\n',SNRvolmlrecov,SNRvol20mlrecov)


%% Generate 1's data
% 2D data vectors
% [8 x numG] = [counts xy(2) nNorm var(4)] x numG  
fprintf('\nBegin image reconstruction ... \n')
fprintf('\n\nGenerate sensitivity/randoms data points ... \n')
tic
% PHI     = -pi/2:pi/sensAngles:pi/2*(1-1/sensAngles);
% indxSup = find(sup1);
% nData   = numel(indxSup);
% nPhi    = numel(PHI);
% data1   = zeros([ 8  nData*nPhi ],'single');
% for kLM = 1:nData
%     r0Tmp = PHI-pi/2;
%     tmp0  = [YXmm(indxSup(kLM)) ; XYmm(indxSup(kLM))];
% 
%     xTmp  = tmp0(1) + sysRad*cos(PHI);
%     yTmp  = tmp0(2) + sysRad*sin(PHI);
% 
%     mTmp = ( yTmp - tmp0(2) ) ./ ( xTmp - tmp0(1) );
%     iTmp = yTmp - mTmp.*xTmp;
%     bTmp = 2*iTmp.*mTmp ./ ( 1 + mTmp.^2 );
%     cTmp = ( iTmp.^2 - sysRad.^2 ) ./ ( 1 + mTmp.^2 );
%     x1Tmp = (-bTmp + sqrt(bTmp.^2 - 4*cTmp) ) / 2;
%     x2Tmp = (-bTmp - sqrt(bTmp.^2 - 4*cTmp) ) / 2;
%     y1Tmp = mTmp.*x1Tmp + iTmp;
%     y2Tmp = mTmp.*x2Tmp + iTmp;
% 
%     mTmpy = ( xTmp - tmp0(1) ) ./ ( yTmp - tmp0(2) );
%     iTmp = xTmp - mTmpy.*yTmp;
%     bTmp = 2*iTmp.*mTmpy ./ ( 1 + mTmpy.^2 );
%     cTmp = ( iTmp.^2 - sysRad.^2 ) ./ ( 1 + mTmpy.^2 );
%     y1Tmpy = (-bTmp + sqrt(bTmp.^2 - 4*cTmp) ) / 2;
%     y2Tmpy = (-bTmp - sqrt(bTmp.^2 - 4*cTmp) ) / 2;
%     x1Tmpy = mTmpy.*y1Tmpy + iTmp;
%     x2Tmpy = mTmpy.*y2Tmpy + iTmp;
% 
%     flipXY        = find(abs(mTmp)>1);
%     x1Tmp(flipXY) = x1Tmpy(flipXY);
%     x2Tmp(flipXY) = x2Tmpy(flipXY);
%     y1Tmp(flipXY) = y1Tmpy(flipXY);
%     y2Tmp(flipXY) = y2Tmpy(flipXY);
% 
%     rTmp = sqrt( (tmp0(1) - x1Tmp).^2 + (tmp0(2) - y1Tmp).^2);
%     lTmp = sqrt( (tmp0(1) - x2Tmp).^2 + (tmp0(2) - y2Tmp).^2);
% 
%     flipLR        = find(rTmp>lTmp);
%     xTmp          = x1Tmp(flipLR);
%     x1Tmp(flipLR) = x2Tmp(flipLR);
%     x2Tmp(flipLR) = xTmp;
%     yTmp          = y1Tmp(flipLR);
%     y1Tmp(flipLR) = y2Tmp(flipLR);
%     y2Tmp(flipLR) = yTmp;
%     dTmp          = rTmp(flipLR);
%     rTmp(flipLR)  = lTmp(flipLR);
%     lTmp(flipLR)  = dTmp;
% 
%     r11Tmp = cos(r0Tmp);
%     r12Tmp = sin(r0Tmp);
%     r21Tmp = -r12Tmp;
%     r22Tmp = r11Tmp;
% 
%     % % Locally invariant (LI) resolution model
%         % Non-collinearity resolution loss
%     k2Tmp = double(lTmp.^2 ./ sin(pi - fwhmAcol/2).^2);
%     bTmp  = double(rTmp.^2 - (lTmp - rTmp).^2 - k2Tmp );
%     cTmp  = double(rTmp.^2 .* (lTmp - rTmp).^2);
%     acolTmp = single(2*(-bTmp - sqrt(bTmp.^2 - 4*cTmp))*sigmaConst^2);
%         % Detector resolution loss
%     vDet1Tmp   = ( sigmaDetmm * (1 - rTmp./lTmp) ).^2;
%     vDet2Tmp   = ( sigmaDetmm *      rTmp./lTmp  ).^2;
%         % Combine variances
%     varDetTmp  = acolTmp + vDet1Tmp + vDet2Tmp;
% 
%     iV11Tmp = r11Tmp.*r11Tmp./varDetTmp + r12Tmp.*r12Tmp./(sigmaTmm^2);
%     iV12Tmp = r21Tmp.*r11Tmp./varDetTmp + r22Tmp.*r12Tmp./(sigmaTmm^2);
%     iV21Tmp = iV12Tmp;
%     iV22Tmp = r21Tmp.*r21Tmp./varDetTmp + r22Tmp.*r22Tmp./(sigmaTmm^2);
% 
%     varNormTmp = 1./(2*pi*sqrt(varDetTmp.*sigmaTmm^2));
% 
%     kStart = (kLM-1)*nPhi+1;
%     kEnd = kStart+nPhi-1;
%     data1(1,kStart:kEnd)  = 1; 
%     data1(2,kStart:kEnd)  = tmp0(1);
%     data1(3,kStart:kEnd)  = tmp0(2);    
%     data1(4,kStart:kEnd)  = varNormTmp;
%     data1(5,kStart:kEnd)  = iV11Tmp;
%     data1(6,kStart:kEnd)  = iV12Tmp;
%     data1(7,kStart:kEnd)  = iV21Tmp;
%     data1(8,kStart:kEnd)  = iV22Tmp;
% 
%     % Detector points
%     x1sens(:,kStart:kEnd) = [x1Tmp; y1Tmp];
%     x2sens(:,kStart:kEnd) = [x2Tmp; y2Tmp];
% 
% end
toc

%% Add attenuation correction
% fprintf('\nGet attenuation for sensitivity image ... \n')
% 
% tic
% sensImAtt = getAttenuationMEX([x1sens(1,:)' -x1sens(2,:)'],[x2sens(1,:)' -x2sens(2,:)'],single(PETFOV),single(PETFOV),single(3.75),flipud(single(muCT)),int32(matSize),int32(matSize),int32(length(x1sens(1,:))));
% toc
% 
% % FLIPPING TO ORIENT TO ROSS'S WAY
% x1sens              = flipud(x1sens);
% x1sens(2,:)         = x1sens(2,:);
% x2sens              = flipud(x2sens); 
% x2sens(2,:)         = x2sens(2,:);
% data1(2:3,:)        = flipud(data1(2:3,:));
% data1(3,:)          = data1(3,:);


%% LI Box AVX2 mex test
tic
% fprintf('\nGenerate sensitivity image and initialize image (Box+Chunk+AVX2)  ... \n')
% % % % % IS THIS HOW WE ADD ATTENUATION??? YES I THINK SO
% data1(1,:) = data1(1,:) .* sensImAtt';
% 
% AT1 = backProjLInv_BoxChunk(data1, ...
%     single(xyGridmm),single(yxGridmm),single(mm2vox),boxStride,int32(32));

fprintf('\nLoad sensitivity image and initialize image (Box+Chunk+AVX2)  ... \n')
load('report_results\AT1_phantom_r125_sysfov300.mat')
AT1 = supDisk.*AT1;
AT1tmp = AT1;
AT1tmp(AT1tmp<1E-5)=1;
f = numG*sup1/nSup1;
data = dataVec;
toc

% save("report_results\AT1_phantom_r125_sysfov300.mat", 'AT1')


%% MLEM
fprintf('\nBegin MLEM (Box+Chunk+AVX2) ... \n')
Af = dataVec; 
for iRecon = 1:numReconIts
    % Positron range
    f = imPosRangeFilt_tripleImage2Dfaster(f, padSize,...
        lungKernel, softKernel, boneKernel, maskLungSoft, maskSoftBone,...
        mixFactorLungSoft, mixFactorSoftBone);
    % Forward project step: An*f+gamma
    Aftmp = forwardProjLInv_BoxChunk(dataVec,f, ...
        single(xyGridmm),single(yxGridmm),single(mm2vox),boxStride,int32(32));
    % Add attenuation
    Aftmp = Aftmp .* muVec';
    % Calcualte ratio: 1/(An*f+gamma)
    g2Af = 1./(Aftmp+1E-4);
    Af(1,:) = g2Af;

    % Add attenuation
    Af(1,:) = Af(1,:) .* muVec';
    % Backproject ratio ATn(1/(An*f+gamma))
    ATg2Af = backProjLInv_BoxChunk(Af, ...
        single(xyGridmm),single(yxGridmm),single(mm2vox),boxStride,int32(32));
    % Compute preconditioner
    preA = f./AT1tmp;
    % Perfrom update: f = (f/AT1)*AnT(1/(An*F+gamma))
    f = preA.*ATg2Af;
    f = f .*sup1;
    % Positron range
    f = imPosRangeFilt_tripleImage2Dfaster(f, padSize,...
        lungKernel, softKernel, boneKernel, maskLungSoft, maskSoftBone,...
        mixFactorLungSoft, mixFactorSoftBone);
    fprintf('%3d ',iRecon)
    if (mod(iRecon,10)==0), fprintf(' ... \n'), end
    figure, imshow(f,[])
    title(sprintf('Iter %d | CTR: %.0f ps | Perp. FWHM: %.1f.', iRecon, (sysCTR.*1e12), fwhmDetPE));

    % saving image for later
    ImRecons{iRecon,kRes} = f;
end
fprintf('\tComplete MLEM (Box+Chunk+AVX2)  ... ')
toc


















