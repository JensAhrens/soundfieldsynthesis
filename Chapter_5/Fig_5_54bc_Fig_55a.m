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

figure_number = '5.54(bc)'; % creates Fig. 5.54(b) and (c)
%figure_number = '5.55(a)'; % creates Fig. 5.55(a)

if ( strcmp( figure_number, '5.55(a)' ) )
    t = 0.013;
else
    t = 0;
end

v     = 150; 
f     = 500; 
c     = 343; 
M     = v/c;

omega = 2*pi*f;

alpha_rot = pi/6;

x_lim = [ -3 3 ];

if ( strcmp( figure_number, '5.55(a)' ) )
    y_lim = [ 0 6 ];
else
    y_lim = [ -3 3 ];
end

% create spatial grid
X        = linspace( x_lim(1), x_lim(2), 1000 );
Y        = linspace( y_lim(1), y_lim(2), 1000 );
[ x, y ] = meshgrid( X, Y );
z        = 0;

% Eq. (5.58)
Delta = sqrt( ( x - v*t ).^2 + ( y.^2 + z.^2 ) .* ( 1 - M^2 ) );

% Eq. (5.57)
tau = ( M .* ( x - v*t ) + Delta ) ./ ( c * ( 1 - M^2 ) );

t_tilde = t - tau;

x_s_tilde = v * t_tilde;

% Eq. (5.75a)
alpha_tilde = atan2( y, x - x_s_tilde );

%%%% convolution with sign (Eq. (5.79)) is ommitted for simplicity %%%% 
    
% Eq. (5.79)
s = cos( alpha_tilde - alpha_rot ) .* 1 / (4*pi) .* exp( 1i .* omega .* t_tilde ) ./ Delta;
                          
figure;

imagesc( X, Y, real( s ), [ -.1 .1 ] );
turn_imagesc;
axis square;
colormap gray;

hold on;
plot( [ x_lim( 1 ) x_lim( 2 ) ], [ 0 0 ], 'k:', 'LineWidth', 2 ); % for dipole field only
hold off;

xlabel( 'x (m)' )
ylabel( 'y (m)' )

graph_defaults;

if ( strcmp( figure_number, '5.54(bc)' ) )

    figure;
    imagesc( X, Y, 20*log10( abs( s ) ), [ -50 -10 ] );
    turn_imagesc;
    axis square;
    colormap gray;
    revert_colormap;
    colorbar;

    hold on;
    plot( [ x_lim(1) x_lim(2) ], [ 0 0 ], 'k:', 'LineWidth', 2 );
    hold off;

    xlabel( 'x (m)' )
    ylabel( 'y (m)' )

    graph_defaults;

end
