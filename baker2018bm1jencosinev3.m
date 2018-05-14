% Nicholas F. Baker
% Benchmark Competition 1
%    Jensen Cosine implementation
% Begun 7 Mar 2018
% Completed 12 Mar 2018 (Original assignment)
% Modified 09 Apr 2018 (For Me575 group project)
% Modified 11 May 2018 (Bug fixes and Cos/Partial Wake code separation)

close all, clear all

nNumCycles = 100;
% NREL 5MW parameters
nNumRtrs = 9; % Number of rotors in field
D = 126.4;                   % Rotor Diameter, needed for constraints

RtrLoc = zeros(nNumRtrs,2);    % Of format (xLocation,yLocation) for each rotor
OptFarm = zeros(nNumCycles,(nNumRtrs*2));
AEPlist = zeros(nNumCycles,1);
BestAEP = [1,0];                   % To store the location of our best AEP. (index, value)

%x0 =  [-500., -500., -500., 0., -500., 500., 0., -500., 0., 0., 0., 500., 500., -500.,500., 0., 500., 500.];
%%{
% Run it and save our results
for i = 1:nNumCycles
    i                    % Output which iteration we're on
    %rng('shuffle');      % Ensures we're using "random" numbers
    rng('default');     % For repeatable "random" numbers
    x0 = ((rand(1,18)-0.5)*1000); % Move between -1000 and 1000
    [xopt,fopt,~,~] = optimize_BenchmarkJens(x0, nNumRtrs);
    
    OptFarm(i,:) = xopt;
    AEPlist(i) = fopt;
    if (BestAEP(2) > AEPlist(i))    % If our new AEP is better (Remember negative switches)
        BestAEP(2) = AEPlist(i);    % Save it
        BestAEP(1) = i;             % And the index of which run we're on
    end
end
%}    
% Plots the optimized rotor locations
for i = 1:nNumRtrs
    %RtrLoc(i,1) = x0((i*2)-1);   % x-coordinate
    %RtrLoc(i,2) = x0(i*2);       % y-coordinate
    RtrLoc(i,1) = OptFarm(BestAEP(1),((i*2)-1));   % x-coordinate
    RtrLoc(i,2) = OptFarm(BestAEP(1),(i*2));       % y-coordinate
end
%RtrLoc

format long
fprintf('%.8f ', BestAEP(2))
%plot(RtrLoc(:,1),RtrLoc(:,2), 'kO', 'MarkerFaceColor', 'k', 'MarkerSize', 20)
plot(RtrLoc(:,1),RtrLoc(:,2), 'kO', 'MarkerSize', 20)
axis([-10*D 10*D -10*D 10*D])
hold off

function [xopt, fopt, exitflag, output] = optimize_BenchmarkJens(x0, nNumRtrs)

    %rng('shuffle');      % Ensures we're using "random" numbers
    rng('default');     % For repeatable "random" numbers
    
    % ------------Starting point and bounds------------
    %{
    x0 =  [-500., -500., ...
           -500., 0., ...
           -500., 500., ...
           0., -500., ...
           0., 0., ...
           0., 500., ...
           500., -500., ...
           500., 0., ...
           500., 500.];
     %}
    %x0 = ((rand(18,1)-0.5).*1000).'; % Move between -1 and 1
       % Of format:
       %  [x1, y1, ...
       %   x2, y2, ... 
    ub = [];
    lb = [];

    % ------------Linear constraints------------
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    %{
        % For plotting wind strength by location
        D = 126.4;  % Rotor Diameter
        Vmin = 7; % (m/s)
        Vmax = 11.4; % (m/s)
        Sigma = 4*D;    % Standard Dev = 4 rotor diameters
        [xgrid,ygrid] = meshgrid(-10*D:1:10*D,-10*D:1:10*D);

        yowza = Vmin + (Vmax - Vmin)* exp( -0.5.*((xgrid./Sigma).^2 + (ygrid./Sigma).^2) );
        figure(2)

        [C,h] = contour(xgrid,ygrid,yowza,6:0.5:13,'blue');
        clabel(C,h,'Labelspacing',250);
        hold on
    %}
    % ------------Objective and Non-linear Constraints------------
    function [f, c, ceq] = objcon(x)
        %--- Analysis Variables ---%
        % AEP variables
        DaysOp = 365;               % Annual days of operation
        HrsOp = 24;                 % Daily hours of operation
        AnlHrsOp = DaysOp * HrsOp;  % Annual hours of operation
        % Wake variables
        Alpha = 0.1;                 % From (Jensen,83)
        % NREL 5MW parameters
        D = 126.4;                   % Rotor Diameter, needed for constraints
        % For degree buckets
        nDegBucket = 1;      % How many degrees to have in each calculated bucket
        nLastDeg = 360 - nDegBucket;
        % Vars for Frequency calc
        degDir = (0:nDegBucket:nLastDeg)';       % A list of every direction to calculate Frequency at
        nNumBuckets = length(degDir);
        FreqWindDir = zeros(nNumBuckets,1);
        [nNumPairs,~] = size(combnk(1:nNumRtrs,2));   % Calculates number of unique pairs for spacing constraints using combinatorics 
        
        %--- Design Variables ---%
        RtrLoc = zeros(nNumRtrs,2);    % Of format (xLocation,yLocation) for each rotor
        for i = 1:nNumRtrs
            RtrLoc(i,1) = x((i*2)-1);   % x-coordinate
            RtrLoc(i,2) = x(i*2);       % y-coordinate
        end
        
        %-- Analysis Functions --%
        % Wind seen from turbine location
        UinfLoc = BenchmarkWindSpeed(RtrLoc(:,1), RtrLoc(:,2));
        % Calculation adjusting wind due to wake
        [uObs] = BenchmarkWakeEffectsFull(UinfLoc, RtrLoc, nDegBucket, Alpha);
        %[uObs] = BenchmarkWakeEffectsPartial(UinfLoc, RtrLoc, nDegBucket, Alpha);
        
        % Power Curve, found in the benchmark paper
        [Pu] = BenchmarkPower(uObs);
        [nNumDir, ~] = size(Pu);
        PwrDir = zeros(nNumDir, 1);
        for i = 1:nNumDir
            PwrDir(i) = sum(Pu(i,:));
        end
        
        % Frequency formula, found in benchmark paper
        nCntr = 0;
        for i = 0:nDegBucket:nLastDeg
            nCntr = nCntr + 1;
            FreqWindDir(nCntr) = BenchmarkWindDist(i);
            %FreqWindDir(nCntr) = integral(@BenchmarkWindDist, i, (i + nDegBucket));
        end      %CheckRose = sum(FreqWindDir);
        
        %{
        % To plot the Wind Frequency distribution
        thetaPlot = deg2rad(0:nDegBucket:(360-nDegBucket));
        polarplot(thetaPlot,FreqWindDir)
        ax = gca;
        ax.ThetaZeroLocation = 'top';
        ax.ThetaDir = 'clockwise';
        %}
        
        %-- Objective Function --%
        % AEP function
        f = -(sum(PwrDir.*FreqWindDir)*AnlHrsOp);  % Multiply the power seen in each direction by the frequency, then sum all power to get total farm power
        
        % Inequality Constraints
        nConstIndex = 0;            % Initialize our index array
        c = zeros(nNumPairs,1);     % Initialize constraint array, for how many unique pairs we have
        for i = 1:(nNumRtrs-1)
            for j = (i+1):nNumRtrs
                nConstIndex = nConstIndex+1;      % Map the constraints to be sequential, in order of pairs
                c(nConstIndex) = (-(RtrLoc(i,1) - RtrLoc(j,1))^2 - (RtrLoc(i,2) - RtrLoc(j,2))^2) + ((2*D)^2);    % Proximity constraint
            end
        end
        %c
        % Equality Constraints (None)
        ceq = [];
    end

    % ------------Call fmincon------------
    options = optimoptions(@fmincon, 'display', 'iter-detailed');
    [xopt, fopt, exitflag, output] = fmincon(@obj, x0, A, b, Aeq, beq, lb, ub, @con, options);
    
    % ------------Separate obj/con (do not change)------------
    function [f] = obj(x)
            [f, ~, ~] = objcon(x);
    end
    function [c, ceq] = con(x)
            [~, c, ceq] = objcon(x);
    end
end
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
function [Pu] = BenchmarkPower(uObs)
    % NREL 5MW parameters
    Prat = 5;       % (MW), Turbine Power rating
    Vci = 3;        % (m/s) Cut-In Wind Speed
    Vrws = 11.4;    % (m/s) Rated Wind Speed
    
    [nNumDir, nNumRotors] = size(uObs);
    Pu = zeros(nNumDir, nNumRotors);
    
    for i = 1:nNumDir
        for j = 1:nNumRotors
            if (uObs(i,j) < Vci)
                Pu(i,j) = 0;
            elseif (Vci <= uObs(i,j)) && (uObs(i,j) <= Vrws)
                Pu(i,j) = Prat * ((uObs(i,j) - Vci) / (Vrws - Vci));
            else    %(uObs > Vrws)
                Pu(i,j) = Prat;
            end
        end
    end
end
function [DistValue] = BenchmarkWindDist(degDir)
    % Given a wind diection (degDir), returns the frequency, as given in
    % Benchmark paper
    
    % Vars from Frequency calc
    Mu(1) = 170;       % (deg) Specified in benchmark paper.
    Mu(2) = 250;       % (deg) Specified in benchmark paper.
    Sig(1) = 20;      % (deg) Specified in benchmark paper.
    Sig(2) = 40;      % (deg) Specified in benchmark paper.
    W(1) = 0.5;    % coefficient to weight frequency. Specified in benchmark paper.
    W(2) = 0.5;    % coefficient to weight frequency. Specified in benchmark paper.
    
    FreqWindDir1stTerm = W(1)*( sqrt(1/ (2*pi*(Sig(1)^2)) ) ) * exp( -( ((degDir-Mu(1)).^2)./(2.*Sig(1)^2) ) );
    FreqWindDir2ndTerm = W(2)*( sqrt(1/ (2*pi*(Sig(2)^2)) ) ) * exp( -( ((degDir-Mu(2)).^2)./(2.*Sig(2)^2) ) );
    DistValue = FreqWindDir1stTerm + FreqWindDir2ndTerm;
end
function [RotatedRtrLoc] = RotatePoints(RtrLoc, Phi)
    % Takes the rotors in RtrLoc and rotates them according to the
    % direction passed as Phi (in degrees, not radians)
    
    Xs = RtrLoc(:,1);
    Ys = RtrLoc(:,2);
    [PolarRtrLoc(:,1),PolarRtrLoc(:,2)] = cart2pol(Xs,Ys);                                      % Translate to polar coordinates
    PolarRtrLoc(:,1) = PolarRtrLoc(:,1) + (deg2rad(90) + deg2rad(Phi));                         % Shift all locations by given direction so wind comes from -x direction (270 deg)
    [RotatedRtrLoc(:,1),RotatedRtrLoc(:,2)] = pol2cart(PolarRtrLoc(:,1), PolarRtrLoc(:,2));     % Translate back to cartesian coord system
end
function [PartialArea, PartialCenter] = JensenPartialWakeOverlap(PriRtrXy, TgtRtrXy, Diam, Alpha)
    %-- Takes two rotated rotor coordinates and determines the area of wake overlap --%
    d = abs(TgtRtrXy(2) - PriRtrXy(2));         % Subtract the y-coordinates for the 'd' distance.
    xDiff = abs(TgtRtrXy(1) - PriRtrXy(1));     % Subrtact the x-coordinates for the x-distance.
    r = Diam/2;                                 % Radius of rotor = half diameter of NREL 5MW rotors
    R = r + (Alpha * xDiff);                    % Radius of wake after having spread to the TgtRtr's location
    
    %- Area calculations -% 
    FirstTerm = (r^2) * acos((d^2 + r^2 - R^2) / (2*d*r));
    SecondTerm = (R^2) * acos((d^2 + R^2 - r^2) / (2*d*R));
    Radical = (1/2) * sqrt( (-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R) );
    PartialArea = FirstTerm + SecondTerm + Radical;
    
    %- Center calculations -%
    PartYcFromPriYc = ((d^2 - r^2 + R^2)/ (2*d));   % Y-coord calculation in lab book
    PartialCenter = zeros(2,1);                     % Initialize our point array. Format: (x,y).
    PartialCenter(1) = TgtRtrXy(1);                 % Partial's x-coord lies in plane of TgtRtr
    if (PriRtrXy(2) < TgtRtrXy(2))                  % If TgtRtr is above the PriRtr
        PartialCenter(2) = PriRtrXy(2) + PartYcFromPriYc;    % Add the distance to center
    else                                            % If TgtRtr is below
        PartialCenter(2) = PriRtrXy(2) - PartYcFromPriYc;    % Subtract it
    end
end
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
function [OrderList] = TurbineOrder(RtrLoc)
    % Determines turbine order from left to right along the x axis.
    % If two turbines have same x-coord then turbine with lower y-coord is
    % ordered first
    nNumRtrs = length(RtrLoc);              % Get the number of turbines
    
    RtrLoc = [RtrLoc, (1:nNumRtrs)'];       % Pair it with a tag for thier position
    OrderedMatrix = sortrows(RtrLoc, 1);  	% Sort according to x position
    OrderList = OrderedMatrix(:,3);         % Get just the order tags, and return it.
end
function [ULoc] = BenchmarkWindSpeed(x,y)
    % Given an (x,y) location, returns the wind speed that rotor sees
    % Equation given in Benchmark Paper (eq 2)
    % Farm Properties
    D = 126.4;  % Rotor Diameter
    Vmin = 7; % (m/s)
    Vmax = 11.4; % (m/s)
    Sigma = 4*D;    % Standard Dev = 4 rotor diameters
    
   ULoc = Vmin + (Vmax - Vmin)* exp( -0.5.*((x./Sigma).^2 + (y./Sigma).^2) );
end
