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

eta = 1; % for Fig. 5.20(a)
%eta = 4; % for Fig. 5.20(b)

L = 1.4;

% create spatial grid
X        = linspace( -L/2, L/2, eta+2 );
Z        = [ -L/4 L/4 ];
[ x, z ] = meshgrid( X, Z );

amplitudes = ones( size( x ) );

% alternate amplitudes of vibrating sections 
for l = 1 : 2 : eta+2
    amplitudes( :, l ) = -1;   
end


figure;

colormap( [ 0 0 0; 1 1 1 ] );

hold on;
% plot plate
h = surf( x, zeros( size( x ) ), z, amplitudes, 'LineStyle', '-' );
    
% plot x axis
line( [ -1 -L/2 ], [ 0  0 ], [ 0 0 ], 'LineStyle', '-', 'Color', [ .5 .5 .5 ], 'LineWidth', 1 );
line( [ L/2 1 ]  , [ 0  0 ], [ 0 0 ], 'LineStyle', '-', 'Color', [ .5 .5 .5 ], 'LineWidth', 1 );
line( 1, 0, 0, 'Marker', '.', 'Color', [ .5 .5 .5 ], 'LineWidth', 1 );

% plot y axis
line( [ 0 0 ],  [   -1   0 ], [ 0 0 ], 'LineStyle', '-', 'Color', [ .5 .5 .5 ], 'LineWidth', 1 );
line( [ 0 0 ],  [ L/8+.1 1 ], [ 0 0 ], 'LineStyle', '-', 'Color', [ .5 .5 .5 ], 'LineWidth', 1 );
line( 0, 1, 0, 'Marker', '.', 'Color', [ .5 .5 .5 ], 'LineWidth', 1 );

% plot z axis
line( [ 0 0 ], [ 0  0 ], [ -1  -L/4 ], 'LineStyle', '-', 'Color', [ .5 .5 .5 ], 'LineWidth', 1 );
line( [ 0 0 ], [ 0  0 ], [ L/4   1  ], 'LineStyle', '-', 'Color', [ .5 .5 .5 ], 'LineWidth', 1 );
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
