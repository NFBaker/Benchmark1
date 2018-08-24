function [uObs] = BenchmarkWakeEffectsFull(UinfLoc, RtrLoc, nDegBucket, Alpha)
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
        nBucketCtr = nBucketCtr + 1;                    % Increase our counter
        RotatedRtrLoc = RotatePoints(RtrLoc, degCur);   % Rotate our points so we're looking from the wind's frame of reference
        IndxOrder = TurbineOrder(RotatedRtrLoc);        % Get index list of which turbines are leftmost to rightmost.
        WakeEfctFctr = ones(nNumRtrs,1);                % Running total of effect of wakes from other turbines, initialize to 100% (1). Reset for every direction.
        
        % Calculate wake effects from each rotor
        for i = 1:(nNumRtrs-1)                          % Don't check last rotor, it won't effect itself, and nothing is behind it
            for j = (i+1):nNumRtrs                      % Check effect on every other further-downstream rotor
                % Get y-coordinates for Primary rotor's (i) left and right pt
                PriXc = RotatedRtrLoc(IndxOrder(i),1);
                PriYc = RotatedRtrLoc(IndxOrder(i),2);
                
                % Store TgtRtr's center point y-coord in a local variable
                TgtXc = RotatedRtrLoc(IndxOrder(j),1);
                TgtYc = RotatedRtrLoc(IndxOrder(j),2);
                
                % Calculate data for Cosine curve (20 degree window)
                Opp = abs(TgtYc - PriYc);
                Adj = TgtXc - PriXc;                % Difference in x positions, value used to evaluate Jensen
                ThetaJensen = atand(Opp/Adj);       % Finds where on the cosine curve we are

                %-- Wake check --%
                % If Center TgtRtr is in wake of PrimRtr, bInWake = true
                bInWake = ( abs(ThetaJensen) <= 20 );
                % Now that we know where it is, do our wake calculations:
                if bInWake              % If we're in the wake:
                    WakeEfctFctr(IndxOrder(j)) = CalculateWakeEffectCosine(RotatedRtrLoc(IndxOrder(i),:), RotatedRtrLoc(IndxOrder(j),:), WakeEfctFctr(IndxOrder(j)), D, Alpha);
                end  % If we're not in the wake, do nothing (no else statement)
            end
        end
        
        uObs(nBucketCtr,:) = WakeEfctFctr';     % Save the wake effect for all rotors in this direction
    end
    
    uObs = uObs .* (UinfLoc)';                  % Multiply wake effect by U infinity seen at each rotor
end