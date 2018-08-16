%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This work is supplementary material for the book                        %
%                                                                         %
% Jens Ahrens, Analytic Methods of Sound Field Synthesis, Springer-Verlag %
% Berlin Heidelberg, 2012, http://dx.doi.org/10.1007/978-3-642-25743-8    %
%                                                                         %
% It has been downloaded from http://soundfieldsynthesis.org and is       %
% licensed under a Creative Commons Attribution-NonCommercial-ShareAlike  %
% 3.0 Unported License. Please cite the book appropriately if you use     %
% these materials in your own work.                                       %
%                                                                         %
% (c) 2012 by Jens Ahrens                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

y_ref    =  1;
delta_x  = .2;
theta_pw = pi/4;
phi_pw   = pi/2;
fs       = 44100;
c        = 343;
tap      = 1675; % time instant to plot in samples

% create frequency axis
no_of_freq_bins = 2049;
F               = linspace( 0, fs, no_of_freq_bins+1 );
F               = F( 1 : end/2 ).';
F(1)            = F(2); % to avoid numerical instabilities
K               = 2*pi*F / c;

X_0             = -10 : delta_x : 10;
[ x_0, k ]      = meshgrid( X_0, K );
y_0             = 0;
   
% Eq. (A.3a), (A.3b)
k_pw_x = k .* cos( theta_pw ) .* sin( phi_pw );
k_pw_y = k .* sin( theta_pw ) .* sin( phi_pw );

% Eq. (3.79)
D = ( 4*i .* exp( -1i .* k_pw_y .* y_ref ) ) ./ ...
            besselh( 0, 2, k_pw_y .* y_ref ) .* ...
                    exp( -1i .* k_pw_x .* x_0) ;

d = real( ifft( [ D; conj( flipud( D( 2 : end-1, : ) ) ) ] ) );

% move t = 0 to the center of the buffer
d = circshift( d, [ (no_of_freq_bins-1)/2 0 ] );
        
% put zeros around to have some headroom
d = [ zeros( 512, length( X_0 ) ); d; zeros( 512, length( X_0 ) ) ];
  
% normalize
d = d ./ max( abs( d( : ) ) );

% create spatial grid
resolution = 400;
X          = linspace( -2, 2, resolution );
Y          = linspace( -1, 3, resolution );
[ x, y ]   = meshgrid( X, Y );

% initialize s
s = zeros( size( x ) );

% avoid confusion
clear x_0;

% this loop evaluates Eq. (5.69)
for index = 1 : length( X_0 )
    
    x_0 = X_0( index );
    
    r = sqrt( ( x - x_0 ).^2 + ( y - y_0 ).^2 );
    
    % from Eq. (5.69)
    t = ( tap/fs - r./c ) .* fs + 1; % in samples

    % Interpolate the impulse responses to find the values at instances
    % t, which correspond to the spatial locations that we are 
    % interested in. 
    d_reshaped = reshape( interp1( ( 1 : size( d, 1 ) ), d( :, index ), t, 'linear'), ...
                                       resolution, resolution );
                 
    % from Eq. (5.69)
    s( find( r > 0 ) ) = s( find( r > 0 ) ) + ...
        d_reshaped( find( r > 0 ) ) ./ r( find( r > 0 ) );
    
end

figure;
imagesc( X, Y, 20*log10( abs( s ) ), [ -20 10 ] );

hold on;
% plot secondary sources
plot( X_0, 0, 'kx' );
hold off;

turn_imagesc;
colorbar;
colormap gray;
revert_colormap;
axis square;
xlabel( 'x (m)' );
ylabel( 'y (m)' );
graph_defaults;
