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

L     = 8;
Alpha = linspace( 0, 2*pi, 2*L + 1 ); % sampling points along azimuth
Beta  = acos( legpts( L ) );          % sampling points along colatitude

[ alpha, beta ] = meshgrid( Alpha, Beta );

% notional radius of 1
grid = ones( size( alpha ) );

% convert to cartesian coordinates
[ Xm, Ym, Zm ] = sph2cart( alpha, beta - pi/2, real( grid ) );

figure;
colormap( [ 0 0 0 ] )
hold on;

% plot lines, take care of occlusion manually
mesh( Xm( :, 2 : end-1 ), ...
      Ym( :, 2 : end-1 ), ...
      Zm( :, 2 : end-1 ), 'LineStyle', '-' );

% backside
mesh( Xm( 1 : 2, [ end-1 end 1 2 ] ), ...
      Ym( 1 : 2, [ end-1 end 1 2 ] ), ...
      Zm( 1 : 2, [ end-1 end 1 2 ] ), 'LineStyle', '-'); 

% plot dots
plot3( Xm( :, 2 : end-1 ), ...
       Ym( :, 2 : end-1 ), ...
       Zm( :, 2 : end-1 ), 'k.' );

% backside
plot3( Xm( 1 : 2, [ end-1 end 1 2 ] ), ...
       Ym( 1 : 2, [ end-1 end 1 2 ] ), ...
       Zm( 1 : 2, [ end-1 end 1 2 ] ), 'k.' );

hold off;
axis equal;
axis off;

view( -95, 30 )

graph_defaults;
