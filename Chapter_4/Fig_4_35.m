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

subfigure = 'a';
%subfigure = 'b';
%subfigure = 'c';
%subfigure = 'd';

f     = linspace( 0, 3000, 200 ).';
c     = 343;
k     = 2.*pi.*f./c;
k(1)  = k(2); % to avoid numerical instabilities
y_ref =  1;
y_s   = -1;
x_s   =  0;


%%%%%%%%%%%%%%%%%%%%%%%%%% prepare spatial fft %%%%%%%%%%%%%%%%%%%%%%%%%%%%
spatial_interval = [ -5 5 ];
delta_x          = .01; % sampling interval for spatial fft in meters

X = spatial_interval( 1 ) : delta_x : spatial_interval( 2 );
Y = linspace( -1, 3, 201 );

k_x_s = (2*pi) / delta_x; % spatial sampling frequency

% create k_x
k_x    = linspace( 0, k_x_s/2, ( length( X ) + 1 ) / 2 ); % positive frequencies
k_x(1) = k_x(2); % to avoid numerical instabilities
k_x    = [ -fliplr( k_x( 2 : end ) ), k_x ]; % adds negative frequencies

% create 2D grid
[ k_x_m, k_m ] = meshgrid( k_x, k );
%%%%%%%%%%%%%%%%%%%%%%% end prepare spatial fft %%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize D_kx
D_kx = zeros( size( k_x_m ) );

indices = find( abs( k_x_m ) <= k_m );

% Eq. (4.52), (4.53), first case
D_kx( indices ) = exp( 1i .* k_x_m( indices ) .* x_s ) .* ...
    besselh( 0, 2, sqrt( k_m( indices ).^2 - k_x_m( indices ).^2 ) * ( y_ref - y_s ) ) ./ ...
        besselh( 0, 2, sqrt( k_m( indices ).^2 - k_x_m( indices ).^2 ) .* y_ref );
    
indices = find( abs( k_x_m ) > k_m );

% Eq. (4.52), (4.53), second case
D_kx( indices ) = exp( 1i .* k_x_m( indices ) .* x_s ) .* ...
    besselk( 0, sqrt( k_x_m( indices ).^2 - k_m( indices ).^2 ) * ( y_ref - y_s ) ) ./ ...
       besselk( 0, sqrt( k_x_m( indices ).^2 - k_m( indices ).^2 ) * y_ref );
    
%%%%%%%%%%%%%%%%%% bandwidth limitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( subfigure == 'c' )
    B_low  = -25;
    B_high =  25;
elseif ( subfigure == 'd' )    
    B_low  =  0;
    B_high = 50;
end
    
if ( subfigure == 'c' || subfigure == 'd' )    
    central_bin = round( length( D_kx ) / 2 );

    indices = [ 1 : central_bin + B_low, central_bin + B_high : length( D_kx ) ];
    D_kx( :, indices ) = 0;

end
%%%%%%%%%%%%%%%%%% end bandwidth limitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% discretization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( subfigure == 'b' || subfigure == 'c' || subfigure == 'd' )    

    D = ifftx( D_kx, [], 2 );

    indices = 1 : length( D );
    indices = indices( mod( indices-1, 20 )~=0 );
    
    % use every 20th value; this correponds to a sampling step of 20 cm 
    D( :, indices ) = 0;
    
    D_kx = fftx( D, [], 2 );

    % normalize a bit to make it look nicer ;)
    D_kx = D_kx ./ abs( D_kx( 20, (end+1)/2 + 5 ) ).* 0.1276;

end
%%%%%%%%%%%%%%%%%%%%%%%% end discretization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
imagesc( k_x, f, 20*log10( abs( D_kx ) ) )
turn_imagesc;
caxis( [ -40 20 ] );
xlim(  [ -80 80 ] );
xlabel( 'kx (rad/m)' );
ylabel( 'f (Hz)' );
colormap gray;
revert_colormap;
colorbar;
axis square;
graph_defaults;
