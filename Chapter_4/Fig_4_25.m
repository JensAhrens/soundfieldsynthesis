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

% both Fig. 4.25(a) and (b) will be created

clear;

M_display = 80;
M         = 27; 
R         = 1.5;
L         = 56;
theta_pw  = pi/2;
phi_pw    = pi/2;
    
f    = linspace( 0, 3000, 200 ).';
k    = (2.*pi.*f)./343;
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

    D_ring_mS( :, m+M_display+1 ) = sum( D .* exponent, 2 );
end

% initlialize G_ring_m
G_ring_m = zeros( length(f), 2*M_display+1 );

% this loop evaluates Eq. (2.34)
for m = -M_display : M_display
    
   disp( [ 'Calculating order ' num2str( m ) ' of ' num2str( M_display ) '.' ] );
    
    for n = abs(m) : M_display + 10 % go towards infinity

        G_ring_m( :, m+M_display+1 ) = G_ring_m( :, m+M_display+1 ) + ...
                    (-1i) .* k .* sphbesselh( n, 2, k.*R ) .* sphharm( n, -m, pi/2, 0 ) .* ...
                            sphbesselj( n, k.*R ) .* sphharm( n, m, pi/2, 0 );
        
        % Note that setting m = 0 in the last sphharm function removes the
        % exponential function. This is easier than typing Eq. (2.34)
        % explicitely...
    end
    
end

% Fig. 4.25(a)
figure;

imagesc( -M_display : M_display, f, 20*log10( abs( G_ring_m .* D_ring_mS ) ), [ -40 -10 ] ); 

turn_imagesc;
colormap gray;
revert_colormap;
colorbar;
axis square;
xlabel( 'm'      );
ylabel( 'f (Hz)' );
graph_defaults;

G_ring_m( :, 1 : -M + M_display )      = 5*eps; % avoid log of 0
G_ring_m( :, M + M_display + 2 : end ) = 5*eps; % avoid log of 0

% Fig. 4.25(b)
figure;

imagesc( -M_display : M_display, f, 20*log10( abs( G_ring_m .* D_ring_mS ) ), [ -40 -10 ] ); 

turn_imagesc;
colormap gray;
revert_colormap;
colorbar;
axis square;
xlabel( 'm'      );
ylabel( 'f (Hz)' );
graph_defaults;
