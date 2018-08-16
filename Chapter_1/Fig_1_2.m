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

% amplitude difference in dB
d = 0; % for Fig. 1.2 a
%d = 6; % for Fig. 1.2 b

a = .085; % head radius

f = 1000; % frequency
c = 343;
k = 2*pi*f/c;

% create a spatial grid
resolution = 300;
X          = linspace( -1.5, 1.5, resolution );
Y          = linspace(  -.5, 2.5, resolution );
[ x, y ]   = meshgrid( X, Y );

%%%% positions of loudspeakers %%%
R   = 2; % distance
x_r = R * cos( 60 / 180 * pi );
y_r = R * sin( 60 / 180 * pi );

x_l = -R * cos( 60 / 180 * pi );
y_l =  R * sin( 60 / 180 * pi);

% sound field of right loudspeaker
S = 10^( d / 20 ) .* point_source( x, y, x_r, y_r, 0, k );
% sound field of left loudspeaker
S = S + point_source( x, y, x_l, y_l, 0, k );

% normalize to center of plot
S = S ./ abs( S( end/2, end/2 ) );

% no sound field inside head
r          = sqrt( x.^2 + y.^2 ); 
S( r < a ) = 0;

figure;
imagesc( X, Y, real( S ) );

hold on;
% plot head
x_circ = linspace( -a, a, 100 );
plot( x_circ,  sqrt( a.^2 - x_circ.^2 ), 'k', 'LineWidth', 2 );
plot( x_circ, -sqrt( a.^2 - x_circ.^2 ), 'k', 'LineWidth', 2 );
% mark loudspeaker positions
plot( [ x_r x_l ], [ y_r y_l ], 'kx' );
hold off;

caxis( [ -.2 .2 ] );
colormap gray; 
turn_imagesc;
axis square;

xlabel( 'x (m)' );
ylabel( 'y (m)' );

graph_defaults;


