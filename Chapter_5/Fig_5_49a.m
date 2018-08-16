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

v   = 600;
f_s = 500; 
c   = 343; 
M   = v / c;

x = 1;
y = 1;
z = 0;

% time instant when Mach cone arrives
t_mach = 1/v * ( x + y * sqrt( M^2 - 1 ) );
t = linspace( t_mach, 0.01, 500 );

% Eq. (5.58)
Delta = sqrt( ( x - v*t ).^2  + ( y.^2 + z.^2 ) .* ( 1 - M^2 ) );

% Eq. (5.72)
f_1 = f_s * ( 1 + M / ( 1 - M^2 ) * ( M + ( x - v*t ) ./ Delta ) );
f_2 = f_s * ( 1 + M / ( 1 - M^2 ) * ( M - ( x - v*t ) ./ Delta ) );

figure;
plot( t * 1000, f_1, 'k', 'LineWidth', 2 )
hold on;
plot( t * 1000, f_2, 'Color', [ .7 .7 .7 ], 'LineWidth', 2 )
hold off;

xlim( [ 3 10 ] )
ylim( [ -20000 20000 ] )
axis square;
grid on;

xlabel( 't (ms)' )
ylabel( 'f (Hz)' )
legend( 'f_1', 'f_2' )

graph_defaults;
