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

function [ S ] = transfer_function_circular_array_single_position( x, f, M )
% this function calculates the transfer function of a dicrete circular
% secondary source distribution to the point x;
% for fullband synthesis, the higher frequencies are calculated using WFS
% to avoid numerical instabilities

if ( M > 27 )
    fc    = 2000;
    f_ex  = f( f < fc ); % for exact solution
    f_wfs = f( f > fc ); % for WFS solution
else
    f_ex  = f;
    f_wfs = [];
end

L = 56;

theta_pw = -pi/2;
phi_pw   = pi/2;
R        = 1.5;
c        = 343;
k        = 2.*pi.*f_ex ./ c;

% initialize D_ring
D_ring = zeros( length( f_ex ), 2*M+1 );

% this loop evaluates Eq. (3.47)
for m = -M : M

    % Eq. (2.38)
    S_breve = 4*pi * 1i.^( -abs(m) ) .* sphharm( abs(m), -m, phi_pw, theta_pw );
    
    % Eq. (2.37a)
    G_breve = -1i .* k .* sphbesselh( abs(m), 2, k.*R ) .* sphharm( abs(m), -m, pi/2, 0 );   

    D_ring( :, m + M + 1 ) = 1 / ( 2*pi*R ) .* S_breve ./ G_breve;

end

alpha_0 = linspace( 0, 2*pi, L+1 );  
alpha_0 = alpha_0( 1 : end-1 );
beta_0  = pi/2;

% initialize S
S_ex = zeros( length( f_ex ), 1 );

% loop over secondary sources
for l = 1 : L 
    
    D_ex = zeros( length( f_ex ), 1 );   
    
    x_0 = R * cos( alpha_0( l ) ) * sin( beta_0 );
    y_0 = R * sin( alpha_0( l ) ) * sin( beta_0 );
    r   = sqrt( ( x(1) - x_0 )^2 + ( x(2) - y_0 )^2 );
    
    % this loop evaluates Eq. (3.49)
    for m = -M : M   
        D_ex = D_ex + D_ring( :, m + M + 1 ) .* exp( 1i * m * alpha_0( l ) );                
    end
    
    % Eq. (1.4)
    S_ex = S_ex + D_ex .* 1 / (4*pi) .* exp( -1i .* k .* r ) ./ r;

end

% normalize manually
S_ex = S_ex ./ 5.9414;

if ( ~isempty( f_wfs ) )

    k = 2.*pi.*f_wfs ./ c;
    
    alpha_n = alpha_0 + pi;
    
    S_wfs = zeros( length( f_wfs ), 1 );
    
    % loop over secondary sources
    for l = 1 : L

        x_0 = R .* cos( alpha_0(l) );
        y_0 = R .* sin( alpha_0(l) );
        r   = sqrt( ( x(1) - x_0 )^2 + ( x(2) - y_0 )^2 );

        % Eq. (3.89)
        if ( cos( alpha_n(l) - theta_pw ) < 0 ) 
            continue;
        end

        D_wfs = cos( theta_pw - alpha_n(l) ) .* sqrt( 1i .* k ) .* ...
                     exp( -1i .* k .* ( R * cos( alpha_0(l) - theta_pw ) ) );

        % Eq. (1.4)
        S_wfs = S_wfs + D_wfs .* 1 ./ (4*pi) .* exp( -1i .* k .* r ) ./ r;

    end

    % normalize manually
    S_wfs = S_wfs ./ .972;
    
else
    S_wfs = [];
    
end

S = [ S_ex; S_wfs ];

end
