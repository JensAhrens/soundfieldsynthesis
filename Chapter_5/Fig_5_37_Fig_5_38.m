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

% Note that this file uses the doublefactorial function from the specfun
% package by Paul Godfrey 
% ( http://www.mathworks.com/matlabcentral/fileexchange/978 ). The
% corresponding file doublefactorial.m can be found in Common\misc\ .

clear;

figure_number = '5.37'; % for Fig. 5.37(a) and (b)
%figure_number = '5.38'; % for Fig. 5.38(a) and (b)

N = 10;

f = 1000;
c = 343;
k = 2*pi*f/c;

% create spatial grid
resolution = 200;
Y          = linspace( -2, 2, resolution );
Z          = linspace( -2, 2, resolution );
[ z, y ]   = meshgrid( Z, Y );
x          = zeros( size( y ) );

% convert to spherical coordinates
r     = sqrt( x.^2 + y.^2 + z.^2 );
alpha = atan2( y, x );
beta  = acos( z ./ r );

% initialze S
S = zeros( size( y ) );

% this loop evaluates Eq. (5.47) with (5.46) or (5.51)
for n = 0 : N - 1
    
    if ( strcmp( figure_number, '5.37' ) )
    
        %%%%% Eq. (5.46) %%%%%
        if ( n == 0 )
            integral = -1;
        elseif ( rem( n, 2 ) == 1 )
            integral = - 1i^( n - 1 ) * doublefactorial( n ) / ( n * ( n + 1 ) * doublefactorial( n - 1 ) );
        elseif ( rem( n, 2 ) == 0 )
            integral = 0;
            continue;
        end
        
    else
        %%%%% Eq. (5.51) %%%%%
        if ( n == 0 )
            integral = -1/2;
        elseif ( n == 1 )
            integral = -1/3;
        elseif ( rem( n, 2 ) == 0 )
            integral = - 1i^n * doublefactorial( n - 1 ) / ( ( 2*n + 1 ) * doublefactorial( n - 2 ) ) * ...
                                                ( ( n + 1 ) / ( ( n + 2 ) * n ) - 1 / ( n - 1 ) );
        elseif ( rem( n, 2 ) == 1 )
            integral = 0;
            continue;
        end
        
    end
     
    % from Eq. (5.47)
    S = S + 2*1i*sqrt(4*pi) .* k .* integral .* ... 
            1i^( -n ) .* sphharm( n, 0, beta, alpha ) .* sphbesselj( n, k.*r ) .* sqrt( 2*n + 1 );
end

% normalize to figure center
S = S ./ max( abs( S( end/2, end/2 ) ) ) * 2;

% Fig. 5.37(a) or Fig. 5.38(a), respectively
figure;

imagesc( Z, Y, real( S ), [ -1 1 ] );

hold on;
% plot vertical line
plot( [ 0 0 ], [ -2 2 ], 'k:' )
% mark focused source
plot( 0, 0, 'kx' )
hold off

colormap gray;
turn_imagesc;
axis square;

xlabel( 'z (m)' );
ylabel( 'y (m)' );

graph_defaults;

% Fig. 5.37(b) or Fig. 5.38(b), respectively
figure;

imagesc( Z, Y, 20 * log10( abs( S ) ), [ -30 10 ] );

hold on;
% plot vertical line
plot( [ 0 0 ], [ -2 2 ], 'k:' )
% mark focused source
plot( 0, 0, 'kx' )
hold off

colormap gray;
%revert_colormap;
colorbar;
turn_imagesc;
axis square;

xlabel( 'z (m)' );
ylabel( 'y (m)' );

graph_defaults;
