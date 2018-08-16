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

f = linspace( 10, 20000, 2000 ).' ;

figure;

S = transfer_function_linear_array_single_position( [ 0 10 ], f );
semilogx( f, 20*log10( abs( S ) ), 'Color', [ 0 0 0 ] );

hold on;

S = transfer_function_linear_array_single_position( [ 0 1 ], f );
semilogx( f, 20*log10( abs( S ) ), 'Color', [ .4 .4 .4 ] );

S = transfer_function_linear_array_single_position( [ 0 .1 ], f );;
semilogx( f, 20*log10( abs( S ) ), 'Color', [ .7 .7 .7 ] );

hold off;

xlim( [ 10 20000 ] );
ylim( [ -30 30 ] );
axis square;
grid on;
xlabel( 'f (Hz)' );
legend( 'y = 10 m', 'y = 1 m', 'y = 0.1 m', 'Location', 'SouthWest' );
graph_defaults;
