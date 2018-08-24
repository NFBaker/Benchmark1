function [RotatedRtrLoc] = RotatePoints(RtrLoc, Phi)
    % Takes the rotors in RtrLoc and rotates them according to the
    % direction passed as Phi (in degrees, not radians)
    
    Xs = RtrLoc(:,1);
    Ys = RtrLoc(:,2);
    [PolarRtrLoc(:,1),PolarRtrLoc(:,2)] = cart2pol(Xs,Ys);                                      % Translate to polar coordinates
    PolarRtrLoc(:,1) = PolarRtrLoc(:,1) + (deg2rad(90) + deg2rad(Phi));                         % Shift all locations by given direction so wind comes from -x direction (270 deg)
    [RotatedRtrLoc(:,1),RotatedRtrLoc(:,2)] = pol2cart(PolarRtrLoc(:,1), PolarRtrLoc(:,2));     % Translate back to cartesian coord system
end