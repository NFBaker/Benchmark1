function [uObs] = BenchmarkWakeEffectsPartial(UinfLoc, RtrLoc, nDegBucket, Alpha)
    %Computes the observed wind at each turbine location due to:
    %       1) Wind strength at location
    %       2) Wake effects per direction
    % In an array matrix where the columns represent every wind direction,
    % as determined by nDegBuckets. uObs is of dimension NumRotorsx(360/nDegBuckets)
    
    Phi = atan(Alpha);                              % Wake angle (radians) is based on this value (Jensen, 83)
    nNumRtrs = length(RtrLoc);                      % Get the number of turbines
    D = 126.4;                                      % Rotor Diameter
    r0 = D/2.;                                      % Rotor radius
    RtrArea = pi * r0^2;                            % Area of the rotor face
    degDir = (0:nDegBucket:(360-nDegBucket))';      % A list of every direction to calculate Frequency at
    nBucketCtr = 0;                                    % Counter for which bucket we're in. Max is 360.
    uObs = zeros(length(degDir), length(RtrLoc));
    
    % Going through each direction from every direction
    for degCur = 0:nDegBucket:(360-nDegBucket)
        nBucketCtr = nBucketCtr + 1;                          % Increase our counter
        RotatedRtrLoc = RotatePoints(RtrLoc, degCur);   % Rotate our points so we're looking from the wind's frame of reference
        IndxOrder = TurbineOrder(RotatedRtrLoc);        % Get index list of which turbines are leftmost to rightmost.
        WakeEfctFctr = ones(nNumRtrs,1);                % Running total of effect of wakes from other turbines, initialize to 100% (1). Reset for every direction.
        
        % Calculate wake effects from each rotor
        for i = 1:(nNumRtrs-1)                          % Don't check last rotor, it won't effect itself, and nothing is behind it
            for j = (i+1):nNumRtrs                      % Check effect on every other furtehr-downstream rotor
                % Get y-coordinates for Primary rotor's (i) left and right pt
                PriXc = RotatedRtrLoc(IndxOrder(i),1);
                PriYc = RotatedRtrLoc(IndxOrder(i),2);
                PriYr = PriYc + r0;
                PriYl = PriYc - r0;
                
                % Store TgtRtr's center point y-coord in a local variable
                TgtXc = RotatedRtrLoc(IndxOrder(j),1);
                TgtYc = RotatedRtrLoc(IndxOrder(j),2);
                TgtYr = TgtYc + r0;
                TgtYl = TgtYc - r0;
                
                % Difference in x positions, value used to evaluate Jensen
                xDiff = abs(TgtXc - PriXc);
                
                % Calculate wake lines from primary at Tgt x val
                PriRlineY = tan(Phi) * (xDiff) + PriYr;   % y - Location of Primary rotor's Right wake limit
                PriLlineY = tan(-Phi) * (xDiff) + PriYl;  % y - Location of Primary rotor's Left wake limit
                
                %-- Partial Wake Check --%
                % If EITHER the left or right points of TgtRtr are in wake of PrimRtr, bInWake = true
                bInWake = ( ((TgtYr <= PriRlineY) ...
                           & (TgtYr >= PriLlineY)) ...
                          | ((TgtYl <= PriRlineY) ...
                           & (TgtYl >= PriLlineY)) );
                % If BOTH the left or right points of TgtRtr are in wake of PrimRtr, bFullWake = true
                bFullWake = ( ((TgtYr <= PriRlineY) ...
                             & (TgtYr >= PriLlineY)) ...
                            & ((TgtYl <= PriRlineY) ...
                             & (TgtYl >= PriLlineY)) );

                % Now that we know where it is, do our wake calculations:
                if bInWake              % If we're in the wake:
                    if bFullWake        % If we're fully waked, do normal jensen calculations
                        WakeEfctFctr(IndxOrder(j)) = CalculateWakeEffectTopHat(RotatedRtrLoc(IndxOrder(i),:), RotatedRtrLoc(IndxOrder(j),:), WakeEfctFctr(IndxOrder(j)), D, Alpha);
                    else    % If we're waked and the Tgt turbine isn't fully in, we're "partially waked", so do the math.
                        % Get the area and center of the overlap slice, referred to as "partial" or "part"
                        [WakedArea,WakedCenter] = JensenPartialWakeOverlap(RotatedRtrLoc(IndxOrder(i),:), RotatedRtrLoc(IndxOrder(j),:), D, Alpha);
                        UnWakedArea = RtrArea - WakedArea;      % Get the area of the rotor remaining
                  %PROBLEM IS IN CALCULATING AREA OF WAKED PORTION
                        PrcntAreaPart = WakedArea/RtrArea;   % Calculate the percent area of the overlap part
                        PrcntAreaMain = UnWakedArea/RtrArea;   % Calculate the percent area of the rest
                        % Calculate weighted wake effects
                        WakeEffectPart = CalculateWakeEffectTopHat(RotatedRtrLoc(IndxOrder(i),:), WakedCenter, 1, D, Alpha);
                        WakeEffectMain = 1; % No wake effect for the main body (100% freestream)
                        WakeEffectWeighted = (PrcntAreaMain * WakeEffectMain) + (PrcntAreaPart * WakeEffectPart);   % The resultant wake effect seen by the rotor, weighted for the percentage in the wake
                        
                        WakeEfctFctr(IndxOrder(j)) = WakeEfctFctr(IndxOrder(j)) * WakeEffectWeighted;
                    end
                end  % If we're not in the wake, do nothing (no else statement)
            end
        end
        
        uObs(nBucketCtr,:) = WakeEfctFctr';     % Save the wake effect for all rotors in this direction
    end
    
    uObs = uObs .* (UinfLoc)';                  % Multiply wake effect by U infinity seen at each rotor
end