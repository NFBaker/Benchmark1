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