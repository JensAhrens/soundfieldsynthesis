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

% Note that for simplicity, this script does not employ Eq. (4.21) but
% rather evaluates Eq. (4.20) numerically to derive D_ring_(m,S).

clear;

subfigure = 'a'; % for Fig. 4.13(a) 
%subfigure = 'b'; % for Fig. 4.13(b)
%subfigure = 'c'; % for Fig. 4.13(c)
%subfigure = 'd'; % for Fig. 4.13(d)

if ( subfigure == 'a' || subfigure == 'b' )
    M = 80; 
else 
    M = 27; 
end

if ( subfigure == 'b' || subfigure == 'd' )
    L = 56; 
end 

M_display = 80; % display range
R         = 1.5;
theta_pw  = pi/2;
phi_pw    = pi/2;

f    = linspace( 0, 3000, 200 ).';
k    = (2.*pi.*f) ./ 343;
k(1) = k(2); % to avoid numerical instabilities

% initialize D_ring_m
D_ring_m = zeros( length(f), 2*M_display+1 );

% this loop evaluates Eq. (3.47)
for m = -M : M
    
    % from Eq. (2.38)
    S_ring_m =  4*pi .* 1i.^( -abs(m) ) .* sphharm( abs(m), -m, phi_pw, theta_pw );

    % from Eq. (2.37a)
    G_ring_m = -1i .* k .* sphbesselh( abs(m), 2, k.*R ) .* sphharm( abs(m), -m, pi/2, 0 );   

    D_ring_m( :, m+M_display+1 ) = 1/(2*pi*R) .* S_ring_m ./ G_ring_m;
    
end

if ( subfigure == 'b' || subfigure == 'd' )
    
    alpha_0 = (2*pi) / L .* ( 0 : L-1 );
    
    % initialize S
    D = zeros( length(f), L );

    % this loop evaluates Eq. (3.48)
    for l = 1 : L

        exponent  = repmat( exp( 1i.* ( -M_display : M_display ) .* alpha_0(l) ), [ length(f) 1 ] );
        D( :, l ) = sum( D_ring_m .* exponent, 2 );

    end

    % from Eq. (4.19)
    D = D./L;

    % initialize D_ring_mS
    D_ring_mS  = zeros( length( f ), 2*M_display+1 );

    % this loop evaluates Eq. (4.20)
    for m = -M_display : M_display
        
        exponent = repmat( exp( -1i .* m .* alpha_0 ), [ length(f) 1 ] );

        D_ring_mS( :, m+M_display+1 ) = sum( D .* exponent,  2);
    end
    
    D_ring_m_abs = 20*log10( abs( D_ring_mS ) );
    
else
    % avoid log of 0
    D_ring_m( D_ring_m == 0 ) = 5*eps;
    
    D_ring_m_abs = 20*log10( abs( D_ring_m ) );
    
end

figure;
imagesc( -M_display : M_display, f, D_ring_m_abs, [ -40 20 ] );

turn_imagesc;
colormap gray; 
revert_colormap;
colorbar;
axis square;
xlabel( 'm'     );
ylabel( 'f (Hz)');
graph_defaults;

