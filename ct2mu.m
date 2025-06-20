% ct2mu.m
function CTAC = ct2mu(CT,CTenergy)


rhoAir       = 0.00128; % g/ml
rhoWater     = 1;       % g/ml
rhoBone      = 1.6;     % g/ml
% Ect          = 70;      % kVp (effective)
[Zair, Rair] = ParseChemicalFormula('Air');
[Zh2o, Rh2o] = ParseChemicalFormula('Water');
Zbone        = [ 1   ;  6   ; 7   ;  8   ; 20   ];     % from literature
Rbone        = [ 3.4 ; 31.4 ; 1.8 ; 36.5 ; 26.8 ]/100; % from literature
 
muAir_CT  = rhoAir   * Rair'  * PhotonAttenuationQ(Zair,CTenergy/1000,'mac')';  % 1/cm
muH2O_CT  = rhoWater * Rh2o'  * PhotonAttenuationQ(Zh2o,CTenergy/1000,'mac')';  % 1/cm
muBone_CT = rhoBone  * Rbone' * PhotonAttenuationQ(Zbone,CTenergy/1000,'mac')'; % 1/cm
 
indxLow   = find(CT <= 0);
indxHigh  = find(CT >  0);
 
muAir_E  = rhoAir   * Rair'  * PhotonAttenuationQ(Zair, 511/1000,'mac')';  % 1/cm
muH2O_E  = rhoWater * Rh2o'  * PhotonAttenuationQ(Zh2o, 511/1000,'mac')';  % 1/cm
muBone_E = rhoBone  * Rbone' * PhotonAttenuationQ(Zbone,511/1000,'mac')';  % 1/cm
   
CTACtmp           = zeros(size(CT));
CTACtmp(indxLow)  = muH2O_E + CT(indxLow)  / 1000 * (muH2O_E  - muAir_E);
CTACtmp(indxHigh) = muH2O_E + CT(indxHigh) / 1000 * (muBone_E - muH2O_E) * ...
    (muH2O_CT - muAir_CT) / (muBone_CT - muH2O_CT);
 
CTACtmp(CTACtmp<muAir_E) = muAir_E;
CTAC                     = CTACtmp/10; % convert to unit [1/mm]

end
    