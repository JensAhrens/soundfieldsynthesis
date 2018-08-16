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

theta_pw = pi/2;
phi_pw   = pi/2;
M        = 79;
R        = 1.5;
d_ref    = 5;

% variables 
f    = linspace( 0, 1000, 200 ).';
f(1) = f(2); % to avoid numerical instabilities
c    = 343;

K = (2*pi*f)/c;

Alpha_0 = linspace( .5*pi, 2.5*pi, 300 );

[ alpha_0, k ] = meshgrid( Alpha_0, K );
%alpha_n     = alpha_0 + pi;
x_0         = R * cos( alpha_0 );
y_0         = R * sin( alpha_0 );


D_exact = zeros( size( k ) );

% this loop evaluates Eq. (3.49) 
for m = -M : M
    disp( [ 'Calculating order ' num2str( m ) ' of ' num2str( M ) '.' ] );
    
    D_exact = D_exact + 1./(2*pi*R) .* ( 4 .* pi .* (-1i).^abs(m) .* sphharm( abs(m), -m, phi_pw, theta_pw ) ) ./ ...
                    ( (-1i) .* k.* sphbesselh( abs(m), 2, k.*R ) .* sphharm( abs(m), -m, pi/2, 0 ) ) .* exp( 1i .* m .* alpha_0 );

end

% Eq. (3.101) reformulated for circular distributions
D_wfs = sqrt( 8 * pi * k * d_ref * 1i ) .* ...
           cos( alpha_0 - theta_pw ) .* ...
                exp( -1i .* k .* cos( alpha_0 - theta_pw ) );

% secondary source selection            
D_wfs( find( alpha_0 < pi   ) ) = 5*eps; % to avoid log of 0                          
D_wfs( find( alpha_0 > 2*pi ) ) = 5*eps; % to avoid log of 0


% Fig. 3.24(a)
figure;
imagesc( Alpha_0 / pi, f, 20*log10( abs( D_exact ) ), [ 0 30 ] );

turn_imagesc;
colorbar;
colormap gray;
revert_colormap;
axis square;
axis tight;
xlabel( '$\alpha$ ($\pi$ rad)', 'Interpreter', 'latex' );
ylabel( 'f (Hz)' );
ylim( [ 0 1000 ] )
graph_defaults;

% Fig. 3.24(b)
figure;
imagesc( Alpha_0 / pi, f, 20*log10( abs( D_wfs ) ), [ 0 30 ] );

turn_imagesc;
colorbar;
colormap gray;
revert_colormap;
axis square;
axis tight;
xlabel( '$\alpha$ ($\pi$ rad)', 'Interpreter', 'latex' );
ylabel( 'f (Hz)' );
ylim( [ 0 1000 ] )
graph_defaults;

