% Nicholas F. Baker
% Benchmark Competition 1
%    Jensen Top Hat/Partial Wake and Cosine/Full Wake implementations
% Begun 7 Mar 2018
% Completed 12 Mar 2018 (Original assignment)
% Modified 09 Apr 2018 (For Me575 group project)
% Modified 11 May 2018 (Bug fixes and Cos/Partial Wake code separation)
% Modified 24 Aug 2018 (Separated repeated functions into their own files)

close all, clear all

global bTopHatOrCos     % Make it so anyone can see
bTopHatOrCos = false;   % True = TopHat and partial wake, False = Cosine and hub wake
%rng('shuffle');        % Ensures we're using "random" numbers
rng('default');         % For repeatable "random" numbers
nNumCycles = 3;
% NREL 5MW parameters
nNumRtrs = 9; % Number of rotors in field
D = 126.4;                      % Rotor Diameter, needed for constraints

RtrLoc = zeros(nNumRtrs,2);     % Of format (xLocation,yLocation) for each rotor
OptFarm = zeros(nNumCycles,(nNumRtrs*2));
AEPlist = zeros(nNumCycles,1);
BestAEP = [1,0];                   % To store the location of our best AEP. (index, value)

%x0 =  [-500., -500., -500., 0., -500., 500., 0., -500., 0., 0., 0., 500., 500., -500.,500., 0., 500., 500.];
%%{
% Run it and save our results
for i = 1:nNumCycles
    i                    % Output which iteration we're on
    x0 = ((rand(1,18)-0.5)*1000); % Move between -1000 and 1000
    [xopt,fopt,~,~] = optimize_BenchmarkJens(x0, nNumRtrs, bTopHatOrCos);
    
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

function [xopt, fopt, exitflag, output] = optimize_BenchmarkJens(x0, nNumRtrs, bTopHatOrCos)

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
        if bTopHatOrCos
            [uObs] = BenchmarkWakeEffectsPartial(UinfLoc, RtrLoc, nDegBucket, Alpha);
        else
            [uObs] = BenchmarkWakeEffectsFull(UinfLoc, RtrLoc, nDegBucket, Alpha);
        end
        
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
        c = ((2*D) - (pdist(RtrLoc, 'euclidean')))';      % Vectorize the distance formula between pairs
        
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