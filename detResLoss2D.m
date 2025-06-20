function [xD1,xD2,xD0,tC] = detResLoss2D(xC1,xC2,xC0,radSim,sigmaDet)

% sigmaConst  = 1/sqrt(8*log(2));
% centIsim    = simFovMat(1)/2; % dont need to do this right?
numT        = numel(xC1(:,1));
% sigmaDet    = sigmaConst*fwhmDet;

% sample detector offsets
offsetDet1  = normrnd(0,sigmaDet,numT,1);
offsetDet2  = normrnd(0,sigmaDet,numT,1);

% sinVec      = (xC1(:,1) - centIsim) ./ radSim;
% cosVec      = (xC1(:,2) - centIsim) ./ radSim;
sinVec      = xC1(:,1) ./ radSim;
cosVec      = xC1(:,2) ./ radSim;
xD1(:,1)    = xC1(:,1) - offsetDet1.*cosVec;
xD1(:,2)    = xC1(:,2) + offsetDet1.*sinVec;

% sinVec      = (xC2(:,1) - centIsim) ./ radSim;
% cosVec      = (xC2(:,2) - centIsim) ./ radSim;
sinVec      = xC2(:,1) ./ radSim;
cosVec      = xC2(:,2) ./ radSim;
xD2(:,1)    = xC2(:,1) - offsetDet2.*cosVec;
xD2(:,2)    = xC2(:,2) + offsetDet2.*sinVec;

tDenomX     = ( xC1(:,1) - xC2(:,1) );
tDenomY     = ( xC1(:,2) - xC2(:,2) );

indxX       = find(abs(tDenomX) >= abs(tDenomY));
indxY       = find(abs(tDenomX) <  abs(tDenomY));
tC          = zeros([numT 1]);
tC(indxX)   = ( xC0(indxX,1) - xC2(indxX,1) ) ./ tDenomX(indxX);
tC(indxY)   = ( xC0(indxY,2) - xC2(indxY,2) ) ./ tDenomY(indxY);

xD0(:,1)    = tC.*xD1(:,1) + (1-tC).*xD2(:,1);
xD0(:,2)    = tC.*xD1(:,2) + (1-tC).*xD2(:,2);

end