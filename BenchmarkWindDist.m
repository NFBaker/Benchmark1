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