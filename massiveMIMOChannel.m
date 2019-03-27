% Federal University of Rio de Janeiro - UFRJ
% Electrical Engineering Program - COPPE
% Signals, Multimedia, and Telecommunications Laboratory - SMT
%
% Author: Rafael da Silva Chaves
% email: rafael.chves@smt.ufrj.br
%
% Abstract: This function generate a multipath channel for multi-user
% multiple-input multiple-output (MU-MIMO) transmission, following a time
% division duplex (TDD) transmission. This channel model is presented in
% [1] and [2].
%
% References:
%
% [1] - A. M. Sayeed, "Deconstructing Multiantenna Fading Channels", IEEE
% Transaction on Signal Processing, vol. 50, no. 10, pp. 2563-2579, Oct.
% 2006
%
% [2] - Y. Gao, W. Jiang and T. Kaiser, "Bidirectional Branch and Bound
% Based Antenna Selection in Massive MIMO Systems" in 2015 IEEE 26th Annual
% International Symposium on Personal, Indoor, and Mobile Radio
% Communications, 2015, pp. 563–568.

function [H, beta, varargout] = massiveMIMOChannel(commcell, ...
                                                   propagation, ...
                                                   fading, ...
                                                   varargin)
% [H, beta, varargout] = multipathMUMIMOChannel(cell, ...
%                                               propagation, ...
%                                               fading, ...
%                                               varargin)
%
% Inputs:
%
%         -- 'array' is a structure with the parameters of the BS antenna
%         array.
%
%           * 'array.type' is the field with the type of the transmitter
%           array, must be a string. Possible values:
%
%               - 'linear'
%
%           * 'array.nAntennas' is the field with the number of antennas in
%           the array, must be a positive integer number.
%
%         -- 'channel' is a strucuture with the parameters of the MU-MIMO
%         channel.
%
%           * 'channel.nUsers' is the field with the number of users in the
%           transmission, must be a positive integer number.
%
%           * 'channel.nPaths' is the field with the number of multipaths
%           for each user, must be a vector with the same size of the
%           number of users.
%
%           * 'channel.angle' is the field with the angles of arrival for
%           each user, must be a cell and each entry of the cell must be a
%           vector with the same size of the multipath number for that
%           user.
%
% Outputs:
%
%         -- 'H' is the uplink channel matrix, it is a complex matrix with
%         size 'array.nAntennas' x 'channel.nUsers'.
%
%         -- 'beta' is the large-scale fading vector, it is a real vector
%         with size 'channel.nUsers'.

% MACROS

N_ARGIN  = 4;                                                              % Number of input arguments
N_ARGOUT = 5;                                                              % Number of output arguments

% Cell parameters

n_antenna             = commcell.nAntennas;                                % Number of transmit antennas at base station
n_user                = commcell.nUsers;                                   % Number of user terminals
R                     = commcell.radius;                                   % Cell's raidus (circumradius) in meters
bs_height             = commcell.bsHeight;                                 % Height of base station in meters
user_height           = commcell.userHeight;                               % Minimum and maximum heights of user terminals in meters
n_path                = commcell.nPath;                                    % Number of multipaths for each user terminal

% Propagation parameters

lambda                = propagation.lambda;                                % Wavelength of transmitted signal
ref_distance          = propagation.refDist;                               % Reference distance for path loss calculation
mean_shadow_fad_dB    = propagation.meanShadowFad;                         % Shadow fading mean in dB
std_dev_shadow_fad_dB = propagation.stdDevShadowFad;                       % Shadow fading standard deviation in dB
path_loss_exponent    = propagation.pathLossExponent;                      % Decay exponent

fading                = upper(fading);                                     % Type of fading that occurs in the transmission

% Testing for errors

% Erros for wrong numbers of input arguments

if (nargin > N_ARGIN)
    error('Wrong number of input arguments');
elseif (nargout > N_ARGOUT)
    error('Wrong number of output arguments');
elseif (nargout > N_ARGOUT-2)
     error('Wrong number of output arguments');
elseif(strcmp(fading,'RICH')     || ...
       strcmp(fading,'RAYLEIGH') || ...
       strcmp(fading,'UR-LOS'))
     error('Invalid fading');
end

% Erros for wrong values in numeric variables

if (n_antenna <= 0)
    error('Number of antennas must be a positive integer number');
elseif (n_user <= 0)
    error('Number of user terminals must be a positive integer number');
elseif (R <=0)
    error('Cell radius must be a positive real number');
elseif (bs_height <= 0)
    error('Base station height must be a positive real number');
elseif (lambda <= 0)
    error('Wavelength must be a positive real number');
elseif (ref_distance <= 0)
    error('Reference distance must be a positive real number');
end

if (size(n_path,1) == 1)
    n_path = repmat(n_path,n_user,1);
elseif (size(n_path,1) ~= n_user)
    error('Invalid size of multipath');
end

r = sqrt(3)/2*R;

H_user = user_height(1) + (user_height(2) - user_height(1))*rand(n_user,1);% Height of user terminals in meters

antenna_spacing = lambda/2;                                                % Antenna spacing of transmitt array

if (nargin == N_ARGIN - 1)
    aux_cord = rand(K,1);
    
    K_1 = sum(aux_cord < 1/3);
    K_2 = sum(aux_cord < 2/3 & aux_cord > 1/3);
    K_3 = sum(aux_cord > 2/3);
    
    u = rand(K,1);
    v = rand(K,1);
    
    u_1 = u(1:K_1,1);
    v_1 = v(1:K_1,1);
    
    u_2 = u(K_1+1:K_1+K_2,1);
    v_2 = v(K_1+1:K_1+K_2,1);
    
    u_3 = u(K_1+K_2+1:K,1);
    v_3 = v(K_1+K_2+1:K,1);

    x_1 = -R/2*u_1 + R*v_1;
    y_1 = r*u_1;
    
    x_2 = -R/2*u_2 - R/2*v_2;
    y_2 = -r*u_2 + r*v_2;
    
    x_3 = R*u_3 - R/2*v_3;
    y_3 = -r*v_3;
    
    x_user = [x_1' x_2' x_3']';
    y_user = [y_1' y_2' y_3']';
elseif (nargin == N_ARGIN)
    coordinate = varargin{1};
    
    x_user   = coordinate.x_user;
    y_user   = coordinate.y_user;
end

theta_user   = atan2(y_user,x_user);                                       % Departure angle in rad

varargout{1} = [x_user y_user];
varargout{2} = theta_user;

% Large-scale Fading Calculation
% The large-scale coefficient is independent of the small-scale effects and
% depende only the distances of BS and users

mean_shadow_fad    = 10^(mean_shadow_fad_dB/20);                           % Shadow fading mean
std_dev_shadow_fad = 10^(std_dev_shadow_fad_dB/20);                        % Shadow fading standard deviation

mu_shadow_fad      = log(mean_shadow_fad^2/sqrt(std_dev_shadow_fad^2 + ...
    mean_shadow_fad^2));
sigma_shadow_fad   = sqrt(log((std_dev_shadow_fad/mean_shadow_fad)^2 + 1));

d_bs_user = sqrt(x_user.^2 + y_user.^2);                                   % Distance between base station and users in meters
r_bs_user = sqrt(d_bs_user.^2 + (bs_height - H_user).^2);                  % Length of the path traveled by the signal in meters

z         = lognrnd(mu_shadow_fad,sigma_shadow_fad,n_user,1);              % Shadow fading
path_loss = ((lambda/(4*pi*ref_distance))^2).* ...
            (ref_distance./r_bs_user).^path_loss_exponent;                 % Path loss
beta      = z.*path_loss;                                                  % Large-scale fading

switch fading
    case 'RICH'
        G = randn(n_antenna,n_user) + 1i*randn(n_antenna,n_user);          % Small-scale fading coefficient matrix
        G = G./sqrt(2);                                                    % Small-scale fading coefficient matrix with unitary variance
        
        H = G*sqrt(diag(beta));                                            % Uplink channel matrix
    case 'RAYLEIGH'
        G = randn(n_antenna,n_user) + 1i*randn(n_antenna,n_user);          % Small-scale fading coefficient matrix
        G = G./sqrt(2);                                                     % Small-scale fading coefficient matrix with unitary variance
        
        H = G*sqrt(diag(beta));                                            % Uplink channel matrix
    case 'UR-LOS'
        steering_vector = steeringVector(n_antenna, ...
                                         theta_user, ...
                                         antenna_spacing, ...
                                         lambda);
        
        g = randn(1,n_user) + 1i*randn(1,n_user);                          % Small-scale fading coefficient vector
        g = g./sqrt(2);                                                     % Small-scale fading coefficient vector with unitary variance
        
        G = repmat(g,n_antenna,1);                                         % Small-scale fading coefficient matrix
        
        H = G.*steering_vector*sqrt(diag(beta));                           % Uplink channel matrix
    case 'SPARSE'
        if (nargin == N_ARGIN - 1)
            x_object = zeros(n_user,n_path(1));
            y_object = zeros(n_user,n_path(1));
                        
            % Generating Interferung object coordinates
            
            for k = 1:n_user
                for n = 1:n_path(k)
                    if(x_user(k) >= 0)
                        x_object(k,n) = x_user(k)*rand(1);
                        y_object(k,n) = -r + 2*r*rand(1);
                        
                        if(x_user(k) > R/2)
                            while(y_object(k,n) + 2*r/R*x_object(k,n) - 2*r > 0 || ...
                                    y_object(k,n) - 2*r/R*x_object(k,n) + 2*r < 0)
                                x_object(k,n) = x_user(k)*rand(1);
                                y_object(k,n) = -r + 2*r*rand(1);
                            end
                        end
                    else
                        x_object(k,n) = x_user(k) - x_user(k)*rand(1);
                        y_object(k,n) = -r + 2*r*rand(1);
                        
                        if(x_user(k) < -R/2)
                            while(y_object(k,n) - 2*r/R*x_object(k,n) - 2*r > 0 || ...
                                    y_object(k,n) + 2*r/R*x_object(k,n) + 2*r < 0)
                                x_object(k,n) = x_user(k) - x_user(k)*rand(1);
                                y_object(k,n) = -r + 2*r*rand(1);
                            end
                        end
                    end
                end
            end
        elseif (nargin == N_ARGIN)
            coordinate = varargin{1};
            
            x_object = coordinate.x_object;
            y_object = coordinate.y_object;
        end
        
        theta_object = atan2(y_object,x_object);                           % Departure angle in rad
        
        varargout{3} = [x_object y_object];
        varargout{4} = theta_object;
        
        G = zeros(n_path(1),n_user);
        H = zeros(n_antenna,n_user);
        
        for k = 1:n_user
            steering_vector = steeringVector(n_antenna, ...
                                             theta_object(k,:)', ...
                                             antenna_spacing, ...
                                             lambda);
            
            G(:,k) = randn(n_path(1),1) + 1i*randn(n_path(1),1); % Small-scale fading coefficient vector
            G(:,k) = G(:,k)./sqrt(2);                                      % Small-scale fading coefficient vector with unitary variance
            
            H(:,k) = sqrt(beta(k))*steering_vector*G(:,k);                 % Uplink channel matrix
        end
    otherwise
end

end