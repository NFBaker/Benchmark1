function [OrderList] = TurbineOrder(RtrLoc)
    % Determines turbine order from left to right along the x axis.
    % If two turbines have same x-coord then turbine with lower y-coord is
    % ordered first
    nNumRtrs = length(RtrLoc);              % Get the number of turbines
    
    RtrLoc = [RtrLoc, (1:nNumRtrs)'];       % Pair it with a tag for thier position
    OrderedMatrix = sortrows(RtrLoc, 1);  	% Sort according to x position
    OrderList = OrderedMatrix(:,3);         % Get just the order tags, and return it.
end