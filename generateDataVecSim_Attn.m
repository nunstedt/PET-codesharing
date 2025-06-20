function [X0mmSim, X1mmSim, X2mmSim, X0rlmmSim, ThetasSim, dataVecSim, muVecSim, idxsNonZeroGT] = ...
    generateDataVecSim_Attn(phantPET_PR, xyGridBoost, yxGridBoost, ...
    muCTSim, PETFOV, matSize, sysRad, numSimAngles)
    
    % chose not to take these as arguments:
    % sigmaConst, fwhmAcol, sigmaDetmm, sigmaTmm 
    
    % Boosted ground truth setup
    idxsNonZeroGT = find(phantPET_PR > 0);                    % Indices of all non-zero voxels in ground truth PET (where activity exists)
    numNonZeroGT = numel(idxsNonZeroGT);                  % Total number of active voxels in PET image    

    X0mmSim        = zeros([2 numNonZeroGT*numSimAngles],'single');
    X1mmSim        = zeros([ 2   numNonZeroGT*numSimAngles ],'single');
    X2mmSim        = zeros([ 2   numNonZeroGT*numSimAngles ],'single');
    % varVecSim      = zeros([ 3   numNonZeroGT*numSimAngles ],'single');
    X0rlmmSim      = zeros([ 2   numNonZeroGT*numSimAngles ],'single');
    ThetasSim      = zeros([ 1   numNonZeroGT*numSimAngles ],'single');
    
    
    % DO WE WANT TO DO IT FOR THE BOOSTED GRID!?!?!?
    % ATM WE DO IT FOR THE BOOSTED GRID!
    [YXmmBoost, XYmmBoost] = meshgrid(xyGridBoost,yxGridBoost); % note x/y swap
    YXmmBoost        = single(YXmmBoost);
    XYmmBoost        = single(XYmmBoost);
    

    fprintf('\nComputing dataVecSim for all non-zero pixels and all angles... \n')
    % STOLEN AND ADAPTED FROM "Generate 1's data"
    tic
    THETAS     = -pi/2:pi/numSimAngles:pi/2*(1-1/numSimAngles);
    dataVecSim   = zeros([ 8  numNonZeroGT*numSimAngles ],'single');
    
    for kLM = 1:numNonZeroGT
        % tic
        r0Tmp = THETAS-pi/2;
        tmp0  = [YXmmBoost(idxsNonZeroGT(kLM)) ; XYmmBoost(idxsNonZeroGT(kLM))];
    
        xTmp  = tmp0(1) + sysRad*cos(THETAS);
        yTmp  = tmp0(2) + sysRad*sin(THETAS);
    
        mTmp = ( yTmp - tmp0(2) ) ./ ( xTmp - tmp0(1) );
        iTmp = yTmp - mTmp.*xTmp;
        bTmp = 2*iTmp.*mTmp ./ ( 1 + mTmp.^2 );
        cTmp = ( iTmp.^2 - sysRad.^2 ) ./ ( 1 + mTmp.^2 );
        x1Tmp = (-bTmp + sqrt(bTmp.^2 - 4*cTmp) ) / 2;
        x2Tmp = (-bTmp - sqrt(bTmp.^2 - 4*cTmp) ) / 2;
        y1Tmp = mTmp.*x1Tmp + iTmp;
        y2Tmp = mTmp.*x2Tmp + iTmp;
    
        mTmpy = ( xTmp - tmp0(1) ) ./ ( yTmp - tmp0(2) );
        iTmp = xTmp - mTmpy.*yTmp;
        bTmp = 2*iTmp.*mTmpy ./ ( 1 + mTmpy.^2 );
        cTmp = ( iTmp.^2 - sysRad.^2 ) ./ ( 1 + mTmpy.^2 );
        y1Tmpy = (-bTmp + sqrt(bTmp.^2 - 4*cTmp) ) / 2;
        y2Tmpy = (-bTmp - sqrt(bTmp.^2 - 4*cTmp) ) / 2;
        x1Tmpy = mTmpy.*y1Tmpy + iTmp;
        x2Tmpy = mTmpy.*y2Tmpy + iTmp;
    
        flipXY        = find(abs(mTmp)>1);
        x1Tmp(flipXY) = x1Tmpy(flipXY);
        x2Tmp(flipXY) = x2Tmpy(flipXY);
        y1Tmp(flipXY) = y1Tmpy(flipXY);
        y2Tmp(flipXY) = y2Tmpy(flipXY);
    
        rTmp = sqrt( (tmp0(1) - x1Tmp).^2 + (tmp0(2) - y1Tmp).^2);
        lTmp = sqrt( (tmp0(1) - x2Tmp).^2 + (tmp0(2) - y2Tmp).^2);
    
        flipLR        = find(rTmp>lTmp);
        xTmp          = x1Tmp(flipLR);
        x1Tmp(flipLR) = x2Tmp(flipLR);
        x2Tmp(flipLR) = xTmp;
        yTmp          = y1Tmp(flipLR);
        y1Tmp(flipLR) = y2Tmp(flipLR);
        y2Tmp(flipLR) = yTmp;
        dTmp          = rTmp(flipLR);
        rTmp(flipLR)  = lTmp(flipLR);
        lTmp(flipLR)  = dTmp;
    
        % r11Tmp = cos(r0Tmp);
        % r12Tmp = sin(r0Tmp);
        % r21Tmp = -r12Tmp;
        % r22Tmp = r11Tmp;
        % 
        % % % Locally invariant (LI) resolution model
        %     % Non-collinearity resolution loss
        % k2Tmp = double(lTmp.^2 ./ sin(pi - fwhmAcol/2).^2);
        % bTmp  = double(rTmp.^2 - (lTmp - rTmp).^2 - k2Tmp );
        % cTmp  = double(rTmp.^2 .* (lTmp - rTmp).^2);
        % acolTmp = single(2*(-bTmp - sqrt(bTmp.^2 - 4*cTmp))*sigmaConst^2);
        %     % Detector resolution loss
        % vDet1Tmp   = ( sigmaDetmm * (1 - rTmp./lTmp) ).^2;
        % vDet2Tmp   = ( sigmaDetmm *      rTmp./lTmp  ).^2;
        %     % Combine variances
        % varDetTmp  = acolTmp + vDet1Tmp + vDet2Tmp;
        % % varDetTmp
        % 
        % iV11Tmp = r11Tmp.*r11Tmp./varDetTmp + r12Tmp.*r12Tmp./(sigmaTmm^2);
        % iV12Tmp = r21Tmp.*r11Tmp./varDetTmp + r22Tmp.*r12Tmp./(sigmaTmm^2);
        % iV21Tmp = iV12Tmp;
        % iV22Tmp = r21Tmp.*r21Tmp./varDetTmp + r22Tmp.*r22Tmp./(sigmaTmm^2);
        % 
        % varNormTmp = 1./(2*pi*sqrt(varDetTmp.*sigmaTmm^2));
        % 


        % indices
        kStart = (kLM-1)*numSimAngles+1;
        kEnd = kStart+numSimAngles-1;
    
        % % Inserting in dataVec
        % % FLIP COORDINATES IS AS IT SHOULD BE HERE!
        dataVecSim(1,kStart:kEnd)  = 1;
        dataVecSim(2,kStart:kEnd)  = tmp0(1);
        dataVecSim(3,kStart:kEnd)  = -tmp0(2);
        % dataVecSim(4,kStart:kEnd)  = varNormTmp;
        % dataVecSim(5,kStart:kEnd)  = iV11Tmp;
        % dataVecSim(6,kStart:kEnd)  = iV12Tmp;
        % dataVecSim(7,kStart:kEnd)  = iV21Tmp;
        % dataVecSim(8,kStart:kEnd)  = iV22Tmp;
        % 
    
        % % Added for attenuation calculation
    
        % FLIP  COORDINATES IS AS IT SHOULD BE HERE!
        X0mmSim(1,kStart:kEnd)  = tmp0(1);
        X0mmSim(2,kStart:kEnd)  = -tmp0(2);
    
        X1mmSim(1,kStart:kEnd)  = x1Tmp;
        X1mmSim(2,kStart:kEnd)  = -y1Tmp;
    
        X2mmSim(1,kStart:kEnd)  = x2Tmp;
        X2mmSim(2,kStart:kEnd)  = -y2Tmp;
   
    
        % r and l
        X0rlmmSim(1,kStart:kEnd)          = rTmp;
        X0rlmmSim(2,kStart:kEnd)          = lTmp;
    
        % Theta
        ThetasSim(:,kStart:kEnd)          = r0Tmp;
    end
    toc
    
    
    fprintf('\nCalculating attenuation for all combinations of pixels and angles in sim data...\n')
    tic
    % Call to MEX function calculation of attenuation (muCT).
    % DONE! TODO: % Check inputs and outputs with erik here - did not find
    % documentation DONE!
    % DONE! TODO 2: % Might want to upsample CT image the same as the PET ground
    % truth... 
    disp("YES")
    muVecSim = getAttenuationMEX( ...
        X1mmSim', ...
        X2mmSim', ...
        single(PETFOV), ...
        single(PETFOV), ...
        single(3.75), ...
        flipud(muCTSim), ...
        int32(size(muCTSim,1)), ...
        int32(size(muCTSim,2)), ...
        int32(length(X1mmSim)));
    
    % muVecSim = single(ones([length(X1mmSim) 1]));
    
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    % load("attenuation_sim\richSim\muVecSim.mat")
    % save("attenuation_sim\richSim\muVecSim.mat", 'muVecSim')
    
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    toc

end