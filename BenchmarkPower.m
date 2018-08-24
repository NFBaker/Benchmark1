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
                Pu(i,j) = Prat * ((uObs(i,j) - Vci) / (Vrws - Vci))^3;
            else    %(uObs > Vrws)
                Pu(i,j) = Prat;
            end
        end
    end
end