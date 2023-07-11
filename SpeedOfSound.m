%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE SPEED OF SOUND %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function a = SpeedOfSound(P,rho)
% Inputs: Pressure and Density

% speed of sound = sqrt(gamma*pressure/rho)
a = sqrt(1.4*P/rho);