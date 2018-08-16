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

window = [ sin( linspace( 0, pi/2, 25 ) ).^2 ...
                            ones( 1, 150 ) ...
                            cos( linspace( 0, pi/2, 25 ) ).^2 ];
                        
figure;
plot( linspace( -2, 2, length( window ) ), window, 'k', 'Color', [ .7 .7 .7 ], 'Linewidth', 2 );
xlabel( 'x (m)' );
ylim( [ 0 1.2 ] );
grid on;
graph_defaults;

