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