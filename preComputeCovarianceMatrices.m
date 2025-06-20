function [dataVec, dataVecG, varMatLImm, iVarMatLImm, varNormLImm, meanVarLImm,varMatGImm] = ...
    preComputeCovarianceMatrices(numG, X0mm, thetas, r, l, ...
    sigmaConst, fwhmAcol, sigmaDetmm, sigmaTmm, ...
    varMatmm, iVarMatmm, normSplatmm)
    disp("YES 1")
    
    tic
    % maxG       = max(phantG(:));
    % X0mm       = zeros([2 numG],'single');
    % X1mm       = zeros([ 2   numG ],'single');
    % X2mm       = zeros([ 2   numG ],'single');
    % X0rlmm     = zeros([ 2   numG ],'single');
    % Theta0     = zeros([ 1   numG ],'single');
    rotMat     = zeros([ 2 2 numG ],'single');
    % % Globally invariant resolution model in mm
    varMatGImm  = zeros([ 2 2 numG ],'single');
    iVarMatGImm = zeros([ 2 2 numG ],'single');
    % % Locally invariant resolution model in mm
    stdVecmm    = zeros([ 2   numG ],'single');
    varMatLImm  = zeros([ 2 2 numG ],'single');
    iVarMatLImm = zeros([ 2 2 numG ],'single');
    varNormLImm = zeros([ 1   numG ],'single');
    
    
    % 2D data vectors
    % [8 x numG] = [counts xy(2) nNorm var(4)] x numG  
    dataVec    = zeros([ 8  numG ],'single');
    dataVecG   = zeros([ 8  numG ],'single');
    meanVarLImm = 0;
    for kLM = 1:numG
    
        tmp0  = X0mm(:,kLM);    
        rTmp = r(kLM);
        lTmp = l(kLM);
    
        r0Tmp = thetas(kLM)-pi/2;
    
        % Rotate 
        tmpRmat         = [ cos(r0Tmp) sin(r0Tmp) ; -sin(r0Tmp) cos(r0Tmp) ];
        tmpW            = 1; %rand(1,'single');  % Count weights, always 1 for data
        rotMat(:,:,kLM) = tmpRmat;
    
        % Globally invariant (GI) resolution model
        varMatGImm(:,:,kLM)  = tmpRmat' * varMatmm  * tmpRmat;
        iVarMatGImm(:,:,kLM) = tmpRmat' * iVarMatmm * tmpRmat;
    
        % Locally invariant (LI) resolution model
            % Non-collinearity resolution loss
        k2Tmp = double(lTmp^2 / sin(pi - fwhmAcol/2)^2);
        bTmp  = double(rTmp^2 - (lTmp - rTmp)^2 - k2Tmp );
        cTmp  = double(rTmp^2 * (lTmp - rTmp)^2);
        acolTmp = single(2*(-bTmp - sqrt(bTmp^2 - 4*cTmp))*sigmaConst^2);
            % Detector resolution loss
        vDet1Tmp   = ( sigmaDetmm * (1 - rTmp/lTmp) )^2;
        vDet2Tmp   = ( sigmaDetmm *      rTmp/lTmp  )^2;
            % Combine variances
        varDetTmp  = acolTmp + vDet1Tmp + vDet2Tmp;
        varMatTmp  = [ varDetTmp        0 ; 0 sigmaTmm^2       ];
        iVarMatTmp = [ 1/varMatTmp(1,1) 0 ; 0 1/varMatTmp(2,2) ];
    
        stdVecmm(:,kLM)      = sqrt([varMatTmp(1,1) varMatTmp(2,2)]);  
        varMatLImm(:,:,kLM)  = tmpRmat' * varMatTmp  * tmpRmat;
        iVarMatLImm(:,:,kLM) = tmpRmat' * iVarMatTmp * tmpRmat;
        varNormLImm(kLM)     = 1/(2*pi*sqrt(det(varMatTmp)));
    
        dataVec(:,kLM)  = [ tmpW ; tmp0 ; varNormLImm(kLM) ; iVarMatLImm(1,1,kLM) ; iVarMatLImm(1,2,kLM); iVarMatLImm(2,1,kLM) ; iVarMatLImm(2,2,kLM) ];
        dataVecG(:,kLM) = [ tmpW ; tmp0 ; normSplatmm ; iVarMatGImm(1,1,kLM) ; iVarMatGImm(1,2,kLM); iVarMatGImm(2,1,kLM) ; iVarMatGImm(2,2,kLM) ];
    
        meanVarLImm = meanVarLImm + varMatTmp/numG;
    
    end
    toc
    
    
    % % if numG<1E5
    % %     figure,
    % %     plot(X1mm(2,:),-X1mm(1,:),'.g'),hold on
    % %     plot(X2mm(2,:),-X2mm(1,:),'.r'),axis square
    % %     xlim(1.*sysRad*[-1 1])
    % %     ylim(1.*sysRad*[-1 1])
    % %     plot(X0mm(2,:),-X0mm(1,:),'.b')
    % %     title("Ross way of plotting")
    % % end
    % % 
    % 
    % % RICHARD PLOTTING STUFF
    % if numG<1E5
    %     figure,
    %     plot(X1mm(1,:),X1mm(2,:),'.g'),hold on
    %     plot(X2mm(1,:),X2mm(2,:),'.r'),axis square
    % 
    %     xlim(1.1*sysRad*[-1 1])
    %     ylim(1.1*sysRad*[-1 1])
    %     axis equal
    %     title("X1,X2 for attenuation check")
    % end

end