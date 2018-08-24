function [WakeEfctFctrFinal] = CalculateWakeEffectTopHat(PriRtrXy, TgtRtrXy, WakeEfctFctrInit, Diam, Alpha)
    %-- Given the current wake effect factor, calculates the new factor
    % given a rotor location --%
    PriXc = PriRtrXy(1);
    TgtXc = TgtRtrXy(1);
    xDiff = abs(TgtXc - PriXc);
    r0 = Diam/2;                        % Rotor radius

    WakeEfctFctrFinal = WakeEfctFctrInit * (1- (2/3)*((r0/(r0 + Alpha*xDiff))^2));  % Multiply by Jensen factor.
end