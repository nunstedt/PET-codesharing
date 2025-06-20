function ind = checkDetectorRadii(X1,X2,X0,sysRad)
    errX1   = vecnorm(X1,2,2)-sysRad;
    errX2   = vecnorm(X2,2,2)-sysRad;
    
    err = [errX1 ; errX2];
    fprintf("\nThe maximum error in detector distance from ring is    %f mm. \n", max(err(:)))
    
    ind = find(vecnorm(X0,2,2) > 0.9*sysRad);
    n   = length(X0(ind,1));
    fprintf("\nThere are %f detections outside or too close to the ring. \n", n);
    
end