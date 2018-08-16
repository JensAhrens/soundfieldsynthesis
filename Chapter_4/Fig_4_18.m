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

subfigure = 'a'; % for Fig. 4.18(a) 
%subfigure = 'b'; % for Fig. 4.18(b)
%subfigure = 'c'; % for Fig. 4.18(c)
%subfigure = 'd'; % for Fig. 4.18(d)

if ( subfigure == 'a' || subfigure == 'c' )
    x_1 = [  0,  0, 0 ]; 
    x_2 = [ .7,  0, 0 ]; 
    x_3 = [  0, .7, 0 ];
    
else
    x_1 = [  0, .69, 0 ]; 
    x_2 = [  0, .7 , 0 ]; 
    x_3 = [  0, .71, 0 ];
    
end

if ( subfigure == 'a' || subfigure == 'b' )
    M = 27;
else
    M = 79; % approx. fullband
end

f = linspace( 10, 20000, 1000 ).';

figure;

S = transfer_function_circular_array_single_position( x_1, f, M );

semilogx( f, 20*log10( abs( S ) ), 'Color', [ 0 0 0 ] );

hold on;

S = transfer_function_circular_array_single_position( x_2, f, M );
    
semilogx( f, 20*log10( abs( S ) ), 'Color', [ .4 .4 .4 ] );

S = transfer_function_circular_array_single_position( x_3, f, M );

semilogx( f, 20*log10( abs( S ) ), 'Color', [ .7 .7 .7 ] );

hold off;

xlim( [  10 20000 ] );
ylim( [ -20    20 ] );
axis square;
grid on;
xlabel( 'f (Hz) ' );

graph_defaults;
