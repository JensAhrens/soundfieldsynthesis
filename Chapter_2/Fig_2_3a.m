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

angle = linspace( 0, 2*pi, 1000 );

figure;

% plot 0th order
out = exp( 1i .* 0 .* angle );
h   = polar( angle, abs( real( out ) ) );
set( h, 'Color', [ 0 0 0 ], 'Linewidth', 2 );

hold on

% plot 1st order
out = exp( 1i .* 1 .* angle );
h   = polar( angle, abs( real( out ) ) );
set( h, 'Color', [ .4 .4 .4 ], 'Linewidth', 2 );

% plot 2nd order
out = exp( 1i .* 2.* angle );
h   = polar(angle, abs( real( out ) ) );
set( h, 'Color', [ .5 .5 .5 ], 'Linewidth', 2 );

% plot 3rd order
out = exp( 1i .* 3 .* angle );
h   = polar(angle, abs( real( out ) ) );
set( h, 'Color', [ .8 .8 .8 ], 'Linewidth', 2 );

axis square;
grid off;
legend( '0', '1', '2', '3', 'Location', 'SouthWest' );
graph_defaults;

