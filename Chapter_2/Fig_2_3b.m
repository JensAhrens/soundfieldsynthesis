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

x = linspace( -1, 1, 1000 ); % argument of legendre functions

figure;

hold on;

out = asslegendre( 1, 0, x );
plot( x, out, 'Color', [ 0 0 0], 'LineWidth', 2 );

out = asslegendre( 2, 1, x );
plot( x, out, 'Color', [ .25 .25 .25 ], 'LineWidth', 2 );

out = asslegendre( 2, 2, x );
plot( x, out, 'Color', [ .4 .4 .4 ], 'LineWidth', 2 );

out = asslegendre( 3, 1, x );
plot( x, out, 'Color', [ .6 .6 .6 ], 'LineWidth', 2 );

out = asslegendre( 3, 2, x );
plot( x, out, 'Color', [ .75 .75 .75 ], 'LineWidth', 2 );

out = asslegendre( 3, 3, x );
plot( x, out, 'Color', [ .9 .9 .9 ], 'LineWidth', 2 );

hold off;

ylim( [ -30 10 ] )
axis square;
xlabel( 'z' );

grid on;

legend( '(1,0)', '(2,1)', '(2,2)', '(3,1)', '(3,2)', '(3,3)', 'Location', 'SouthWest' );

graph_defaults;
