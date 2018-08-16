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

% Note: For convenience, this script evaluates Eq. (4.28) numerically.

clear;

M        = 80;
theta_pw = pi/2;
R        = 1.5;
y_ref    = .75;

% number of secondary sources (integration points), is chosen high enough 
% to avoid considerable spectral repetitions in the region of interest
L        = 400; 
c        = 343;

% set up 2D grid
f       = linspace( 0, 3000, 200 ).';
K       = (2*pi*f)./c;
Alpha_0 = (2*pi)/L .* ( 0 : L-1 );
[ alpha_0, k ] = meshgrid( Alpha_0, K );
x_0     = R * cos( alpha_0 );
y_0     = R * sin( alpha_0 );

% Eq. (4.31)
D = sqrt( 8*pi .* k .* y_ref / i ) .*  cos( theta_pw - alpha_0  ) .* exp( -1i .* k .* cos( theta_pw - alpha_0 ) ) ;

% window w( alpha_0 )
D( find( alpha_0 / pi < 1 ) ) = 0;                          
D( find( alpha_0 / pi > 2 ) ) = 0;

% initialize D_ring_mS
D_ring_mS = zeros( length(f), 2*M+1 );

% this loop evaluates Eq. (4.28)
for m = -M : M
    D_ring_mS( :, m+M+1 ) = sum( D .* exp( -1i .* m .* alpha_0 ), 2 );
end

% from Eq. (4.19)
D_ring_mS = D_ring_mS / L;

% avoid log of 0
D_ring_mS( D_ring_mS == 0 ) = 5*eps;

figure;

imagesc( -M : M, f, 20*log10( abs( D_ring_mS ) ), [ -30 10 ] );
turn_imagesc;
colormap gray; 
revert_colormap;
colorbar;
axis square;
xlabel( 'm'      );
ylabel( 'f (Hz)' );
graph_defaults;
