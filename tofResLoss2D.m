function [deltaT,xF0] = tofResLoss2D(xD1,xD2,tT,sigmaTmm)
%{
deltaT :    assume the first and second photon are from same decay, and
            detected at time t1 and t2, then deltaT are difference between
            t1 and t2
xF0 :       the results after CTR blur
deltaR :    similar with deltaT, the first and second photon are from
            different decays
%}

sigmaConst = 1/sqrt(8*log(2));
cSpeed     = physconst('LightSpeed')*1000;      % Speed of light mm/sec

distD = sqrt( (xD1(:,1)-xD2(:,1)).^2 + (xD1(:,2)-xD2(:,2)).^2 );

numT = numel(tT);
% sample detector offsets
offsetCTR = normrnd(0,sigmaTmm,numT,1) ./ distD;

tPosCTR = tT + offsetCTR;
deltaT = (tPosCTR-0.5).*distD/cSpeed;

xF0(:,1) = tPosCTR.*xD1(:,1) + (1-tPosCTR).*xD2(:,1);
xF0(:,2) = tPosCTR.*xD1(:,2) + (1-tPosCTR).*xD2(:,2);

% distR  = sqrt( (R1(:,1)-R2(:,1)).^2 + (R1(:,2)-R2(:,2)).^2 );
% deltaR = (tR-0.5).*distR/cSpeed;


end
