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

eta = 1; % for Fig. 5.21(a)
%eta = 4; % for Fig. 5.21(b)

a = .6;

% create spatial grid
alpha           = linspace( 0, 2*pi, 2*eta*10 + 1 );
beta            = linspace( 0,   pi, 20 );  
[ alpha, beta ] = meshgrid( alpha, beta );

amplitudes = ones( size( alpha ) );

% alternate amplitudes of vibrating sections 
for l = 1 : 2*10 : 2*eta*10
    amplitudes( :, l:l+10 ) = -1;   
end

[ Xm, Ym, Zm ] = sph2cart( alpha, pi/2-beta, a );

figure;

colormap( [ 0 0 0; 1 1 1 ] )

hold on;
h = surf( Xm, Ym, Zm, amplitudes, 'LineStyle', '-' );
    
% plot x axis
line( [ -1 -a ], [ 0  0 ], [ 0 0 ], 'LineStyle', '-', 'Color', [ .5 .5 .5 ], 'LineWidth', 1 );
line( [  a  1 ], [ 0  0 ], [ 0 0 ], 'LineStyle', '-', 'Color', [ .5 .5 .5 ], 'LineWidth', 1 );
line( 1, 0, 0, 'Marker', '.', 'Color', [ .5 .5 .5 ], 'LineWidth', 1 );

% plot y axis
line( [ 0 0 ], [ -1 -a ], [ 0 0 ], 'LineStyle', '-', 'Color', [ .5 .5 .5 ], 'LineWidth', 1 );
line( [ 0 0 ], [  a  1 ], [ 0 0 ], 'LineStyle', '-', 'Color', [ .5 .5 .5 ], 'LineWidth', 1 );
line( 0, 1, 0, 'Marker', '.', 'Color', [ .5 .5 .5 ], 'LineWidth', 1 );

% plot z axis
line( [ 0 0 ], [ 0 0 ], [ -1 -a ], 'LineStyle', '-', 'Color', [ .5 .5 .5 ], 'LineWidth', 1 );
line( [ 0 0 ], [ 0 0 ], [  a  1 ], 'LineStyle', '-', 'Color', [ .5 .5 .5 ], 'LineWidth', 1 );
line( 0, 0, 1, 'Marker', '.', 'Color', [ .5 .5 .5 ], 'LineWidth', 1 );

hold off;

axis equal;
axis off;

view( 60, 30 )

text( 1.2,  0 ,  0 , '$x$', 'Interpreter', 'latex' );
text(  0 , 1.2,  0 , '$y$', 'Interpreter', 'latex' );
text(  0 ,  0 , 1.2, '$z$', 'Interpreter', 'latex' );

axis( [ -.6 .6 -.6 .6 -.6 .6 ] );

graph_defaults;
