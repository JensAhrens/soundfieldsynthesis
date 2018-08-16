% This script is based on a script by Benny Bernschütz. 

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

alpha_o = pi/6;
beta_o  = pi/2;

% set up spatial grid 
resolution      = 50;
Alpha           = linspace( 0, 2*pi, resolution );  % azimuth
Beta            = linspace( 0,   pi, resolution );  % colatitude
[ alpha, beta ] = meshgrid( Alpha, Beta );

% spherical law of cosines
gamma = acos( cos( alpha - alpha_o ) .* cos( beta - beta_o ) );

% Eq. (5.81)
S_bar = cos( gamma );

% convert colatitude to elevation and then everything to Cartesian
% coordinates
[ Xm, Ym, Zm ] = sph2cart( alpha, pi/2 - beta, abs( real( S_bar ) ) );

figure;

colormap gray;
axes( 'position', [ 0.1 0.2 .8 .68 ] ); 

surf( Xm, Ym, Zm ); 

hold on;
% plot coordinate axes
line( [ -1 1 ], [ .0 .0 ], [  0 0 ], 'Marker', '.', 'LineStyle', '-', 'Color', [ .5 .5 .5 ], 'LineWidth', [ 1 ] );
line( [  0 0 ], [ -1  1 ], [  0 0 ], 'Marker', '.', 'LineStyle', '-', 'Color', [ .5 .5 .5 ], 'LineWidth', [ 1 ] );
line( [  0 0 ], [  0  0 ], [ -1 1 ], 'Marker', '.', 'LineStyle', '-', 'Color', [ .5 .5 .5 ], 'LineWidth', [ 1 ] );
hold off;

axis( [ -1 1 -1 1 -1 1 ] );

view(35, 25);

light( 'Position', [ 1 -5 3 ] ); 
lighting phong; 
camzoom( 1.3 );

xlabel( '$x$', 'Interpreter', 'latex' );
ylabel( '$y$', 'Interpreter', 'latex' );
zlabel( '$z$', 'Interpreter', 'latex' );

graph_defaults;
