function [X1C,X2C,X0C] = nonColl2D(X1,X2,X0,sysRad) % xT1,xT2,xT0,clT,fwhmThetaNonC,simFovMat,radSim
    


    sz = length(X0(:,1));        
    r  = vecnorm(X0-X1,2,2);
    l  = vecnorm(X2-X1,2,2);

    % WHAT SHOULD THE VALUE BE HERE? the value of h matches 0.0022*diameter
    % at center   
    sigma = getNCstd(r,l);
    
    % input the non coll equation!
    h       = normrnd(0,sigma,[sz,1]);
    
    % find 'center of LOR'-vector
    Xc      = 0.5*X1 + 0.5*X2;
    
    % normalize
    XcHat   = Xc ./ vecnorm(Xc,2,2);
    
    % find center vectors of zero length (LOR pass through origin)
    indZero = find(isnan(XcHat(:,1)));
    
    % rotation matrix 90 degrees
    R       = [0 -1 ; 1 0];
    
    % set NaN parts of XcHat orthogonal to LOR
    XcHat(indZero,:) = ((X2(indZero,:)-X1(indZero,:)) ./ vecnorm((X2(indZero,:)-X1(indZero,:)),2,2) );
    for k = 1:length(indZero)
        V = R * XcHat(indZero(k),:)';
        XcHat(indZero(k),:) = V';
        
    end
    
    X0C     = X0 + h .* XcHat;
    
    % get distance to new center and check that it is inside circle rad
    d       = sysRad^2 - vecnorm(Xc+h.*XcHat,2,2).^2;
    Ix      = find( d < 0 );
    X0C(Ix,:) = X0(Ix,:) - h(Ix) .* XcHat(Ix,:);
    d(Ix)   = sysRad^2 - vecnorm(Xc(Ix,:)-h(Ix).*XcHat(Ix,:),2,2).^2;
    d       = realsqrt(d);
    
%     cosVec  = Xc(:,1) ./ vecnorm(Xc,2,2);
%     sinVec  = Xc(:,2) ./ vecnorm(Xc,2,2);
    cosVec  = XcHat(:,1) ./ vecnorm(XcHat,2,2);
    sinVec  = XcHat(:,2) ./ vecnorm(XcHat,2,2);
    
    X1C(:,1) = Xc(:,1) + d.*sinVec;
    X1C(:,2) = Xc(:,2) - d.*cosVec;
    X2C(:,1) = Xc(:,1) - d.*sinVec;
    X2C(:,2) = Xc(:,2) + d.*cosVec;
    
% end
    
    
%     ROSS    
%     sigmaConst  = 1/sqrt(8*log(2));
%     centIsim    = 0 ; % simFovMat(1)/2;
%     theta_nonC  = pi*(180-fwhmThetaNonC)/180;
%     sin2_nonC   = sin(theta_nonC)^2;
%     k2Tmp       = clT.^2 / sin2_nonC;
% 
%     xTc(:,1) = 0.5*(xT1(:,1) + xT2(:,1));
%     xTc(:,2) = 0.5*(xT1(:,2) + xT2(:,2));
%     rTc      = sqrt( (xTc(:,1) - centIsim ).^2 + (xTc(:,2) - centIsim ).^2);
% 
% 
%     % rTmp  = sqrt(radSim^2 - (clT/2).^2);
%     r1Tmp = sqrt( (xT1(:,1)-xT0(:,1)).^2 + (xT1(:,2)-xT0(:,2)).^2 );
%     r2Tmp = sqrt( (xT0(:,1)-xT2(:,1)).^2 + (xT0(:,2)-xT2(:,2)).^2 );
% 
%     fwhmNonC   = sqrt(2*(k2Tmp-r1Tmp.^2-r2Tmp.^2) - 2*sqrt((k2Tmp-r1Tmp.^2-r2Tmp.^2).^2 - 4*r1Tmp.^2.*r2Tmp.^2) );
%     offsetNonC = normrnd(0,sigmaConst*fwhmNonC);
% 
%     cosVec = ( xTc(:,1) - centIsim ) ./ rTc;
%     sinVec = ( xTc(:,2) - centIsim ) ./ rTc;
%     cosVec(isnan(cosVec)) = 1;
%     sinVec(isnan(sinVec)) = 0;
% 
%     xCc(:,1) = centIsim + (rTc + offsetNonC).*cosVec;
%     xCc(:,2) = centIsim + (rTc + offsetNonC).*sinVec;
% 
%     xC0(:,1) = xT0(:,1) + offsetNonC.*cosVec;
%     xC0(:,2) = xT0(:,2) + offsetNonC.*sinVec;
% 
%     lNonC = 2*sqrt( radSim.^2 - (rTc + offsetNonC).^2 );
% 
%     xC1(:,1) = xCc(:,1) - lNonC.*sinVec/2;
%     xC1(:,2) = xCc(:,2) + lNonC.*cosVec/2;
%     xC2(:,1) = xCc(:,1) + lNonC.*sinVec/2;
%     xC2(:,2) = xCc(:,2) - lNonC.*cosVec/2;
    
    
end