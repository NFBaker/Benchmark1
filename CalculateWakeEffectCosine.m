function [WakeEfctFctrFinal] = CalculateWakeEffectCosine(PriRtrXy, TgtRtrXy, WakeEfctFctrInit, Diam, Alpha)
    %-- Given the current wake effect factor, calculates the new factor
    % given a rotor location --%
    PriYc = PriRtrXy(2);
    TgtYc = TgtRtrXy(2);
    PriXc = PriRtrXy(1);
    TgtXc = TgtRtrXy(1);
    xDiff = abs(TgtXc - PriXc);
    r0 = Diam/2;                        % Rotor radius

    Opp = abs(TgtYc - PriYc);
    Adj = TgtXc - PriXc;
    ThetaJensen = atand(Opp/Adj);       % Finds where on the cosine curve we are
    
    fJensen = (1 + cosd(9 * ThetaJensen))/2;     % Taken from (Jensen 83)
    WakeEfctFctrFinal = WakeEfctFctrInit * (1- (2/3)*(((fJensen * r0)/(r0 + Alpha*xDiff))^2));  % Multiply by Jensen factor.
end