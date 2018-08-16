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

figure_number = '5.31(a)'; % for Fig. 5.31(a)
%figure_number = '5.31(b)'; % for Fig. 5.31(b)
%figure_number = '5.32(b)'; % for Fig. 5.32(b)
%figure_number = '5.33(b)'; % for Fig. 5.33(b)

N     = 81;  
f_max = 3000;

f     = linspace( 0, f_max, 200 ).';
f(1)  = f(2); % to avoid numerical instabilities
c     = 343;
omega = 2*pi*f;
k     = omega/c;

alpha_s = pi;
beta_s  = pi/2;
r_s     = 1;

if ( strcmp( figure_number, '5.31(a)' ) )
    r = r_s / 2;
else
    r = 2 * r_s; 
end

% initialize S_ring_n0
S_ring_n0 = zeros( length( k ), N );

for n = 0 : N - 1

    S_ring_n0( :, n+1 ) = (-1i) .* k .* sphbesselh( n, 2, k .* r_s ) .* ...
                            sphharm( n, 0, beta_s, alpha_s ) .* sphbesselj( n, k .* r );  
        
    if ( strcmp( figure_number, '5.32(b)' ) )
        % Eq. (5.39)
        S_ring_n0( f < (n*c)/(2*pi*r_s), n+1 ) = 0; 
    end       
        
end 

if ( strcmp( figure_number, '5.33(b)' ) )

    % this loop applies Eq. (5.40)
    for n = 1 : length( f )
        n_max = floor( ( 2*pi*f( n ) ) / c * r_s ) + 1;
        
        S_ring_n0( n, 1 : n_max ) = S_ring_n0( n, 1 : n_max ) .* ...
                       ( 0.5 * ( cos( ( 1 : n_max ) ./ n_max * pi ) + 1 ) );
        
        S_ring_n0( n, n_max+1 : end ) = 0;
        
    end
    
end

figure;
imagesc( 0 : N - 1, f, 20*log10( abs( S_ring_n0 ) ), [ -60 -20 ] ); 
turn_imagesc;
colormap gray;
revert_colormap;
colorbar;
axis square;
xlabel( 'n' );
ylabel( 'f (Hz)' );

hold on;
% plot n = omega/c r_s boundary
plot( 0 : N - 1, ( 0 : N - 1 ) * c / (2*pi*r_s), 'k', 'Linewidth', 2 );
hold off;

graph_defaults;
