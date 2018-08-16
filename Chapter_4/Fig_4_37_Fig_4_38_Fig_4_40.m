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

subfigure = '4.37(a)';
%subfigure = '4.37(b)';
%subfigure = '4.38(a)';
%subfigure = '4.38(b)';
%subfigure = '4.40(a)';
%subfigure = '4.40(b)';

f     = linspace( 0, 3000, 300 ).';
c     = 343;
k     = 2.*pi.*f./c;
k(1)  = k(2); % to avoid numerical instabilities
y_ref =  1;
y_s   = -1;
x_s   =  0;

theta_pw = 3/8*pi;
phi_pw   = pi/2;

% create 2D grid
k_x            = linspace( -k(end), k(end), 300 ); 
[ k_x_m, k_m ] = meshgrid( k_x, k );

% Eq. (A.3)
k_pw_x = k_m .* cos( theta_pw ) .* sin( phi_pw );
k_pw_y = k_m .* sin( theta_pw ) .* sin( phi_pw );

% Eq. (3.78)
D_kx = ( 4*1i .* exp( -1i .* k_pw_y .* y_ref ) ) ./ besselh( 0, 2, k_pw_y .* y_ref );
                
% mimick the Dirac in Eq. (3.78)
tolerance = .2;                       
D_kx( ( k_x_m < k_pw_x - tolerance ) | ( k_x_m > k_pw_x + tolerance ) ) = 5*eps;

if ( strcmp( subfigure, '4.38(a)' ) || strcmp( subfigure, '4.38(b)' ) || ... 
     strcmp( subfigure, '4.40(a)' ) || strcmp( subfigure, '4.40(b)' ) )

    % bandwidth limitation
    D_kx( abs(k_x_m) > 15 ) = 5*eps;

end

if ( strcmp( subfigure, '4.37(b)' ) || strcmp( subfigure, '4.38(b)' ) || ...
     strcmp( subfigure, '4.40(a)' ) || strcmp( subfigure, '4.40(b)' ) )

    % add repetitions ( Eq. (4.47) )
    delta_x0 = .2;
    period   = 2*pi/delta_x0;

    period_bins = round( period / ( k_x(end) - k_x(1) ) * length( k_x ) );

    D_kx = D_kx + [ D_kx( :, period_bins : end ), zeros( size( D_kx, 1), period_bins-1 ) ];

    D_kx = D_kx + [ D_kx( :, 2*period_bins : end ), zeros( size( D_kx, 1 ), 2*period_bins-1 ) ];

    D_kx = D_kx + [ zeros( size( D_kx, 1 ), period_bins ), D_kx( :, 1 : end-period_bins ) ];

    D_kx = D_kx + [ zeros( size( D_kx, 1 ), 2*period_bins ) D_kx( :, 1 : end-2*period_bins ) ];

end

if ( strcmp( subfigure, '4.40(a)' ) || strcmp( subfigure, '4.40(b)' ) )
 
    % initialize G_kx
    G_kx = zeros( size( k_m ) );

    % Eq. (C.10), first case
    G_kx( abs( k_x_m ) <= k_m ) = -1i/4 * ...
    besselh( 0, 2, sqrt( k_m( abs( k_x_m ) <= k_m ).^2 - k_x_m( abs( k_x_m ) <= k_m ).^2 ) .* y_ref );

    % Eq. (C.10), second case
    G_kx( abs( k_x_m ) > k_m ) = 1/(2*pi) * ...
    besselk( 0, sqrt( k_x_m( abs( k_x_m ) > k_m ).^2 - k_m( abs( k_x_m ) > k_m ).^2 ) .* y_ref );

    if ( strcmp( subfigure, '4.40(b)' ) )
        % Eq. (4.55)
        G_kx( abs(k_x_m) > 15 ) = 5*eps;
    end
    
    % Eq. (3.71)
    S_kx = D_kx .* G_kx;
    
    data_to_plot = S_kx;

else
    
    data_to_plot = D_kx;
    
end

figure;
imagesc( k_x, f, 20*log10( abs( data_to_plot ) ) )
turn_imagesc;

if ( strcmp( subfigure, '4.40(a)' ) || strcmp( subfigure, '4.40(b)' ) )
    caxis( [ -30  0 ] );
else
    caxis( [   0 30 ] );    
end

xlabel( 'kx (rad/m)' );
ylabel( 'f (Hz)' );
colormap gray;
revert_colormap;
colorbar;
axis square;
graph_defaults;
