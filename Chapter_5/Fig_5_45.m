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

subfigure = 'ab'; % for Fig. 5.45(a) and (b)
%subfigure = 'cd'; % for Fig. 5.45(c) and (d)

if ( subfigure == 'ab' )
    v = 120; 
elseif ( subfigure == 'cd' )
    v = 240;
end

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

% Eq. (5.58)
Delta = sqrt( ( x - v*t ).^2  + ( y.^2 + z.^2 ) .* ( 1 - M^2 ) );

% Eq. (5.57)
tau = ( M .* ( x - v*t ) + Delta ) / ( c * ( 1 - M^2 ) );

% Eq. (5.60)
s = 1 / (4*pi) .* exp( 1i .* omega .* ( t - tau ) ) ./ Delta;

% normalize
s = s ./ abs( s( end/2, end/2 ) );

figure;
imagesc( X, Y, real( s ), [ -1.5 1.5 ] );
turn_imagesc;
axis square;
colormap gray;

hold on;
% plot trajectory
plot( [ -3 3 ], [ 0 0 ], 'k:' )
hold off;

xlabel( 'x (m)' )
ylabel( 'y (m)' )

graph_defaults;

figure;
imagesc( X, Y, 20*log10( abs( s ) ) );
caxis( [ -40 20 ] )
turn_imagesc;
axis square;
colormap gray; 
revert_colormap;
colorbar;

hold on;
% plot trajectory
plot( [ -3 3 ], [ 0 0 ], 'k:' )
hold off;

xlabel( 'x (m)' )
ylabel( 'y (m)' )

graph_defaults;
