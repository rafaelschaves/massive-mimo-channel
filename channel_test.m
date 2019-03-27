clear;
close all;
clc;

M = 500;                                                                   % Number of antennas at the base station
K = 5;                                                                     % Number of users at the cell

commcell.nAntennas       = M;                                              % Number of Antennas
commcell.nUsers          = K;                                              % Number of Users
commcell.radius          = 500;                                            % Cell's raidus (circumradius) in meters
commcell.bsHeight        = 30;                                             % Height of base station in meters
commcell.userHeight      = [1 2];                                          % Height of user terminals in meters ([min max])
commcell.nPaths          = 1;                                              % Number of Multipaths
commcell.frequency       = 2e9;                                            % Carrier frequency in Hz 
commcell.meanShadowFad   = 0;                                              % Shadow fading mean in dB
commcell.stdDevShadowFad = 8;                                              % Shadow fading standard deviation in dB
commcell.city            = 'large';                                        % Type of city

[H_ra, beta_ra] = massiveMIMOChannel(commcell,'rayleigh');
[H_ur, beta_ur] = massiveMIMOChannel(commcell,'ur-los');
[H_sp, beta_sp] = massiveMIMOChannel(commcell,'sparse');
    
