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
% (c) 2015 by Jens Ahrens                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

% This script shows how to simulate a moving source with a trajectory that
% is not parallel to the x-axis. The trick is to rotate the underlying 
% spatial grid, do the simulation as before, and then rotate the coordinate 
% system. This example is not contained in the book.
% See http://www.soundfieldsynthesis.org/moving-sound-sources.

v     = 240; % source velocity
theta = 30 / 180 * pi; % rotation angle of source trajectory
f     = 500; 
t     = 0;
c     = 343; 
omega = 2*pi*f;
M     = v / c;

% create spatial grid
X        = linspace( -3, 3, 500 );
Y        = linspace( -1, 5, 500 );
[ x, y ] = meshgrid( X, Y );
z        = 0;

% rotate spatial grid in opposite direction compared to the source 
% trajectory (Trick!)
x_rot = cos( -theta ) .* x - sin( -theta ) .* y;
y_rot = sin( -theta ) .* x + cos( -theta ) .* y;

% Eq. (5.58)
Delta = sqrt( ( x_rot - v*t ).^2  + ( y_rot.^2 + z.^2 ) .* ( 1 - M^2 ) );

% Eq. (5.57)
tau = ( M .* ( x_rot - v*t ) + Delta ) / ( c * ( 1 - M^2 ) );

% Eq. (5.60)
s = 1 / (4*pi) .* exp( 1i .* omega .* ( t - tau ) ) ./ Delta;

% normalize
s = s ./ abs( s( end/2, end/2 ) );

figure;
% plotting with imagesc() inherently rotates spatial grid and in turn 
% the coordinate system back
imagesc( X, Y, real( s ), [ -1.5 1.5 ] );

set( gca, 'YDir','normal' );
axis square;
colormap jet;

hold on;
% plot trajectory
plot( [ -cos( theta ) cos( theta ) ] * 5, [ -sin( theta ) sin( theta ) ] * 5, 'Color', [ 1 1 1 ], 'Linestyle', '--', 'Linewidth', 2 );
hold off;

xlabel( '$x$ (m)', 'Interpreter', 'latex' );
ylabel( '$y$ (m)', 'Interpreter', 'latex' );

set( gcf, 'Color', [ 1, 1, 1 ] );
