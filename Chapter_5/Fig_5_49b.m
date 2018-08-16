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

v     = 600;
f     = 500; 
c     = 343; 
omega = 2*pi*f;
M     = v / c;

x = 1;
y = 1;
z = 0;

% time instant when Mach cone arrives
t_mach = 1/v * ( x + y * sqrt( M^2 - 1 ) );
t      = linspace( t_mach, 0.01, 400 );

% Eq. (5.58)
Delta = sqrt( ( x - v*t ).^2  + ( y.^2 + z.^2 ) .* ( 1 - M^2 ) );

% Eq. (5.64)
tau_1 = ( M .* ( x - v*t ) + Delta ) / ( c * ( 1 - M^2 ) );
tau_2 = ( M .* ( x - v*t ) - Delta ) / ( c * ( 1 - M^2 ) );

% Eq. (5.63) and (5.62)
s = 1 / (4*pi) .* exp( 1i .* omega .* ( t - tau_1 ) ) ./ Delta + ...
    1 / (4*pi) .* exp( 1i .* omega .* ( t - tau_2 ) ) ./ Delta;

figure;
plot( t * 1000, abs( s ), 'k', 'LineWidth', 2' )

xlim( [ 3 10 ] )
ylim( [ 0  2 ] )
axis square;
grid on;

xlabel( 't (ms)' )
ylabel( '| s( x, t ) |' )

graph_defaults;
