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

% Note that the numerical precision of Matlab is not high enough to
% evaluate Eq. (4.14). D_(n,S)^m is therefore derived by evaluating 
% Eq. (4.12) discretely.

clear;

subfigure = '4.7(a)';
%subfigure = '4.7(b)';
%subfigure = '4.9(a)';
%subfigure = '4.9(b)';

if ( strcmp( subfigure, '4.7(a)' ) || strcmp( subfigure, '4.7(b)' ) )
    M = 79;
    
else
    M = 28;
    
end

L         = 28;
N_display = 80; % display range
M_display = 80; % display range 
R         = 1.5;
theta_pw  = -pi/2;
phi_pw    = pi/2;
f_max     = 3000;

f    = linspace( 0, f_max, 200 ).';
k    = (2.*pi.*f)./343;
k(1) = k(2); % to avoid numerical instabilities

% spatial sampling grid
alpha_0             = linspace( 0, 2*pi, 2*L+1 );
alpha_0             = alpha_0( 1:end-1 );
[ beta_0, weights ] = legpts( L );
beta_0              = acos( beta_0 );
weights             = weights.';

% initialize coefficients of continuous driving function
D_ring_nm_cont = zeros( length(k), 2*M_display + 1, N_display );

% this loop calculates Eq. (3.20)
for n = 0 : M
  
    % Eq. (2.37a)
    G_breve = -1i .* k .* sphbesselh( n, 2, k.*R ) .* sphharm( n, 0, 0, 0 );
    
    for m = -n : n 
         S_breve = 4*pi .* (-1i).^n .* sphharm( n, -m, phi_pw, theta_pw );
         
         D_ring_nm_cont( :, m+M_display+1, n+1 ) = ...
                       1/(2*pi*R^2) * sqrt( ( 2*n+1 ) / ( 4*pi ) ) .* S_breve./G_breve;
    end
end

% if sampled driving function is considered
if ( strcmp( subfigure, '4.7(b)' ) || strcmp( subfigure, '4.9(b)' ) )

    % initialize D
    D = zeros( length( k ), 2*( L )^2 );

    index_sec_source = 1;

    fprintf( '\nEvaluating Eq. (3.21) ...\n\n' );

    % this loop evaluates Eq. (3.21)
    for l_1 = 1 : 2*L % loop over all azimuths
        disp( [ 'Calculating azimuth number ' num2str( l_1 ) ' of ' num2str( 2*L ) '.' ] );

        for l_2 = 1 : L % loop over all colatitutes

            for n = 0 : M                               
                for m = -n : n                                
                    D( :, index_sec_source ) = D( :, index_sec_source ) +  ...
                                    D_ring_nm_cont( :, m+M_display+1, n+1 ) .* ...
                                        sphharm( n, m, beta_0( l_2 ), alpha_0( l_1 ) );                
                end
            end

            index_sec_source = index_sec_source + 1;

        end % l_2 
    end % l_1

    % free some memory
    clear D_ring_nm_cont x y z;

    fprintf( '\nEvaluating Eq. (4.12) ...\n\n' );

    % initialize coefficients of sampled driving function
    D_ring_nm_samp = zeros( length( k ), 2*M_display + 1, N_display );

    % this loop evaluates Eq. (4.12) discretely
    for n = 0 : N_display
        disp( [ 'Calculating order ' num2str( n ) ' of ' num2str( N_display ) '.' ] );

        for m = -n : n

            index_sec_source = 1;

            % loop over loudspeakers
            for l_1 = 1 : 2*L % loop over all azimuths
                for l_2 = 1 : L % loop over all colatitudes

                    D_ring_nm_samp ( :, m+M_display+1, n+1 ) = ...
                        D_ring_nm_samp ( :, m+M_display+1, n+1 ) + ...
                            weights( l_2 ) .* D( :, index_sec_source ) .* ...
                                sphharm( n, -m, beta_0( l_2 ), alpha_0( l_1 ) );               

                    index_sec_source = index_sec_source + 1;

                end % l_2
            end % l_1

        end % m
    end % n
    
    % avoid log of 0
    D_ring_nm_samp( D_ring_nm_samp == 0 ) = 5*eps;
    
    D_ring_nm_abs = 20*log10( abs( D_ring_nm_samp ) );
        
    % free some memory
    clear D D_ring_nm_samp;

else
    
    % avoid log of 0
    D_ring_nm_cont( D_ring_nm_cont == 0 ) = 5*eps;
    
    D_ring_nm_abs = 20*log10( abs( D_ring_nm_cont ) );
    
    % free some memory
    clear D_ring_nm_cont;

end


% preparations for plot
n             = 0 : N_display - 1;
m             = -M_display : M_display;
D_ring_nm_abs = shiftdim( D_ring_nm_abs , 1 );


% clip data to make plot nicer
data_min = -30;
data_max =  30; 

D_ring_nm_abs( D_ring_nm_abs < data_min ) = data_min;
D_ring_nm_abs( D_ring_nm_abs > data_max ) = data_max;


f_max = f(end);
n = 0 : N_display-1;
m = -M_display : M_display;

figure;

h = slice( n, m, f, D_ring_nm_abs, [], [], f );

alpha( 'color' );
set( h, 'EdgeColor', 'none', 'FaceAlpha', 'flat' );

% set axes properties
axis normal;
daspect( [ .025 .05 1 ] );
view( -45, 15 );
axis( [ n(1) n(end)+1 m(1) m(end) f(1) f(end) ] );

set( gca, 'XTick', [   0  20 40 60 80 ] );
set( gca, 'YTick', [ -80 -40  0 40 80 ] );
set( gca, 'YDir',  'reverse' );

colormap gray;
revert_colormap;
caxis( [ data_min data_max ] );
colorbar;

box on;

xlabel( '$n$',       'FontSize', 13, 'Interpreter', 'latex' );
ylabel( '$m$',       'FontSize', 13, 'Interpreter', 'latex' );
zlabel( '$f$~(Hz)',  'FontSize', 13, 'Interpreter', 'latex' );

graph_defaults;
