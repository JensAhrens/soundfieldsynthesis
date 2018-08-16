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

% This script shows how to simulate WFS of a moving source with a 
% trajectory that is not parallel to the x-axis. The trick is to rotate the 
% underlying spatial grid as well as the secondary source distribution, do 
% the simulation, and then rotate the coordinate system. The result is a 
% secondary source distribution that is parallel to the x-axis but the 
% source trajectory is not. This example is not contained in the book.
% See http://www.soundfieldsynthesis.org/moving-sound-sources.

v = 240;

theta = 30 / 180 * pi; % rotation angle of source trajectory
f     = 500; 
t     = 0;
c     = 343; 
omega = 2*pi*f;
M     = v / c;
dx_0  = 0.1;
d_ref = 1;

% create spatial grid
X        = linspace( -3, 3, 500 );
Y        = linspace( -1, 5, 500 );
[ x, y ] = meshgrid( X, Y );
z        = 0;

% rotate spatial grid in opposite direction compared to the source 
% trajectory (Trick!)
x_rot = cos( -theta ) .* x - sin( -theta ) .* y;
y_rot = sin( -theta ) .* x + cos( -theta ) .* y;

% initialize s
s = zeros( size( x_rot ) );

x_0_prime = -5 : dx_0 : 5; % this moves along the secondary source distribution
y_0_prime = ones( size( x_0_prime ) );

% rotate secondary source distribution (the source still moves along x-axis)
x_0s = cos( -theta ) * x_0_prime - sin( -theta ) * y_0_prime; 
y_0s = sin( -theta ) * x_0_prime + cos( -theta ) * y_0_prime; 
z_0 = 0;

% loop over secondary sources (Eq. (5.69))
for index = 1 : length( x_0_prime ) 

    x_0 = x_0s( index );
    y_0 = y_0s( index );
    
    disp( [ 'Processing secondary source at x_0 = ' num2str( x_0, '%0.1f' ) ' m.' ] ); 
    
    r_0 = sqrt( ( x_rot - x_0 ).^2 + ( y_rot - y_0 ).^2 + ( z - z_0 ).^2  );
    
    % Eq. (5.58) and (5.69)
    Delta = sqrt( ( x_0 - v * ( t - r_0/c ) ).^2  + ( y_0.^2 + z_0.^2 ) .* ( 1 - M^2 ) );

    % Eq. (5.57) and (5.69)
    tau = ( M .* ( x_0 - v * ( t - r_0/c ) ) + Delta ) / ( c * ( 1 - M^2 ) );

    % Eq. (5.65)
    dx = ( ( x_0 - v*t ) ./ Delta.^2 + 1 / ( c * ( 1 - M^2 ) ) * ( M + ( x_0 - v*t ) ./ Delta ) * 1i*omega ) .* exp( 1i .* omega .* ( t - r_0/c - tau ) ) ./ Delta;
    
    % Eq. (5.66)
    dy =  y_0 ./ Delta .* ( ( 1 - M^2 ) ./ Delta + 1/c * 1i*omega ) .* exp( 1i .* omega .* ( t - r_0/c - tau ) ) ./ Delta;
    
    % Eq. (3.93)
    d = sqrt( 2*pi*d_ref / ( 1i * omega/c ) ) .* 2 .* ( cos( -theta + pi/2 ) * dx + sin( -theta + pi/2 ) * dy );
    
    % Eq. (5.69)
    s = s + 1 / (4*pi) .* d ./ r_0;

end
    
% normalize
s = s ./ abs( s( end/2, end/2 ) );

% plot inherently rotates sound field by -theta (now, the source moves
% along a trajectory that is not parallel to the secondary source
% distribution
figure;
imagesc( X, Y, real( s ), [ -1.5 1.5 ] );

set( gca, 'YDir','normal' );
axis square;
colormap jet;

hold on;
% plot secondary sources
plot( ( -3 : dx_0 : 3 ), 1, 'kx' )
% plot source trajectory
plot( [ -cos( theta ) cos( theta ) ] * 5, [ -sin( theta ) sin( theta ) ] * 5, 'Color', [ 1 1 1 ], 'Linestyle', '--', 'Linewidth', 2 )
% plot position of virtual source
plot( 0, 0, 'k.' )
hold off;

xlabel( '$x$ (m)', 'Interpreter', 'latex' );
ylabel( '$y$ (m)', 'Interpreter', 'latex' );

set( gcf, 'Color', [ 1, 1, 1 ] );


