function [dataVecSampled, X0mmSampled, X1mmSampled, X2mmSampled, X0rlmmSampled, ThetasSampled, muVecSampled, numG] = ...
    sampleDataVecSim_Attn(dataVecSim, X0mmSim, X1mmSim, X2mmSim, ...
    X0rlmmSim, ThetasSim, muVecSim, phantPET_PR, idxsNonZeroGT, gNum, sysRad, numSimAngles)
    disp("YES reportresults")
    fprintf('\nSampling indices from dataVecSim... ')
    tic
    nonZeroPhantPETList = phantPET_PR(idxsNonZeroGT);
    nonZeroPhantPETListtmp = nonZeroPhantPETList(:,ones(1,numSimAngles)); % extending to include angles for each pixel
    nonZeroPhantPETListSim = reshape(nonZeroPhantPETListtmp',[numel(nonZeroPhantPETListtmp) 1]); % flatten for indices to fit loop.
    nonZeroPhantPETAttListSim = muVecSim .* nonZeroPhantPETListSim; 
    
    PDFListSim = gNum .* nonZeroPhantPETAttListSim / sum(nonZeroPhantPETAttListSim(:));
    
    % LOOK OVER NAMES ON STUFF
    rng(42)
    phantGList   = poissrnd(PDFListSim);
    while (sum(phantGList(:))==0)
       phantGList   = poissrnd(PDFListSim);
    end
    toc
    
    fprintf('\nBuilding dataVec... \n')
    tic
    % PLAN:
    % To catch if same pixel and angle (LOR) is sampled twice
    numG             = sum(phantGList(:)); % NUMBER OF SAMPLED EVENTS
    dataVecSampled          = zeros([ 8  numG ],'single'); % .
    X0mmSampled             = zeros([ 2  numG ],'single'); 
    X1mmSampled             = zeros([ 2  numG ],'single'); 
    X2mmSampled             = zeros([ 2  numG ],'single'); 
    muVecSampled            = zeros([ numG  1 ],'single');
    X0rlmmSampled           = zeros([ 2   numG ],'single');
    ThetasSampled           = zeros([ 1   numG ],'single');
    
    maxValPhantGList = max(phantGList); % max number of counts on a LOR
    kStart = 1;
    for numLORCounts = 1:maxValPhantGList
        sampledDataSimIdxsNumLORCounts = find(phantGList >= numLORCounts);
    
        % % dataVec
        tempDataVec = dataVecSim(:,sampledDataSimIdxsNumLORCounts);
        % muVec
        tempMuVec   = muVecSim(sampledDataSimIdxsNumLORCounts,:); % reverse dim order due to output from C-function.
        % X0mm, X1mm, X2mm
        tempX0mm    = X0mmSim(:,sampledDataSimIdxsNumLORCounts);
        tempX1mm    = X1mmSim(:,sampledDataSimIdxsNumLORCounts);
        tempX2mm    = X2mmSim(:,sampledDataSimIdxsNumLORCounts);
        % r and l
        tempX0rlmm = X0rlmmSim(:, sampledDataSimIdxsNumLORCounts);
        % Thetas
        tempThetas = ThetasSim(:,sampledDataSimIdxsNumLORCounts);
    
        % Calc indices
        tempLength = size(tempX0mm, 2);
        kEnd = kStart+tempLength-1;
    
        % Insert
        dataVecSampled(:,kStart:kEnd) = tempDataVec;
        muVecSampled(kStart:kEnd,:)   = tempMuVec;
        X0mmSampled(:,kStart:kEnd)    = tempX0mm;
        X1mmSampled(:,kStart:kEnd)    = tempX1mm;
        X2mmSampled(:,kStart:kEnd)    = tempX2mm;
        X0rlmmSampled(:,kStart:kEnd)  = tempX0rlmm;
        ThetasSampled(:,kStart:kEnd)  = tempThetas;
    
        % Update start index
        kStart = kEnd+1; 
    end
    toc
        
    % % RICHARD PLOTTING STUFF
    % if numG<1E5
    %     figure,
    %     plot(X1mmSampled(1,:),X1mmSampled(2,:),'.g'),hold on
    %     plot(X2mmSampled(1,:),X2mmSampled(2,:),'.r'),axis square
    % 
    %     % plot(dataVec(2,:),dataVec(3,:),'.y') % should be same
    %     plot(X0mmSampled(1,:), X0mmSampled(2,:), '.b')     % as this one
    % 
    %     xlim(1.1*sysRad*[-1 1])
    %     ylim(1.1*sysRad*[-1 1])
    %     axis equal
    %     title("X0,X1,X2 from richSIM")
    % end

end