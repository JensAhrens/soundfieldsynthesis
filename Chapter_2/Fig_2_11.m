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

% refer also to the comment on line 129

subfigure = 'a'; % for Fig. 2.11(a)
%subfigure = 'b'; % for Fig. 2.11(b)
%subfigure = 'c'; % for Fig. 2.11(c)
%subfigure = 'd'; % for Fig. 2.11(d)
%subfigure = 'e'; % for Fig. 2.11(e)
%subfigure = 'f'; % for Fig. 2.11(f)

resolution = 200;

% large scale
if ( subfigure == 'a' )
    % spatial locations to evaluate
    X = linspace( -2, 2, resolution );
    Y = linspace( -2, 2, resolution );
    
    t = 115; % index of time sample to plot 
    
elseif ( subfigure == 'c' )
    % spatial locations to evaluate
    X = linspace( -2, 2, resolution );
    Y = linspace( -2, 2, resolution );
    
    t = 122; % index of time sample to plot 
    
elseif ( subfigure == 'e' )
    % spatial locations to evaluate
    X = linspace( -2, 2, resolution );
    Y = linspace( -2, 2, resolution );
    
    t = 129; % index of time sample to plot 
    
% small scale
elseif ( subfigure == 'b' )
    % spatial locations to evaluate
    X = linspace( -.5, .5, resolution );
    Y = linspace( -.5, .5, resolution );

    t = 115; % index of time sample to plot 
    
elseif ( subfigure == 'd' )
    % spatial locations to evaluate
    X = linspace( -.5, .5, resolution );
    Y = linspace( -.5, .5, resolution );

    t = 122; % index of time sample to plot 
    
elseif ( subfigure == 'f' )
    % spatial locations to evaluate
    X = linspace( -.5, .5, resolution );
    Y = linspace( -.5, .5, resolution );

    t = 129; % index of time sample to plot 
    
end

% order of plane wave
N = 6;

% propagation direction of plane wave
theta_pw = pi/2;
phi_pw   = pi/2;

% create frequency axis 
f_s  = 10000; % sample rate
number_of_frequency_bins = 257; % overall number of sampling points
f    = linspace( 0, f_s, number_of_frequency_bins + 1 ).';
f    = f( 1 : end/2 ); % we calculate only the lower half of the spectrum
f(1) = f(2); % to avoid numerical instabilities

c = 343;
K = 2*pi*f/c;

% create space/frequency grid
[ x, y, k ] = meshgrid( X, Y, K );
r           = sqrt( x.^2 + y.^2 );
kr          = k.*r;
alpha       = atan2( y, x );
beta        = pi/2;

% free some memory
clear x y r k;

% initialize S
S = zeros( length(X), length(Y), length(K) );

% this loop evaluates Eq. (2.38)
for n = 0 : N-1
    disp( [ 'Calculating order ' num2str( n ) ' of ' num2str( N-1 ) '.' ] );
    
    for m = -n : n
        S = S + 4 .* pi .* 1i^( -n ) .* sphharm( n, -m, phi_pw, theta_pw ) .* ...
                    sphbesselj( n, kr ) .* sphharm( n, m, beta, alpha );
        
    end
end

% free some memory
clear kr alpha;

% this will hold the complete spectrum
S_full = zeros( size( S, 1 ), size( S, 2 ), 2*size( S, 3 ) - 1 );

% create upper half of spectrum from lower half
S_full( :, :,           1 : (end+1)/2 ) = S;
S_full( :, :, (end+1)/2+1 :  end      ) = conj( flipdim( S( :, :, 2 : end ), 3 ) );

% free some memory
clear S;

% transfer sound field to time domain
s = ifft( S_full, [], 3 );

% free some memory
clear S_full;

% neglect unwanted imaginary components that might occur due to limited numerical precision
s = real( s );

% logarithmic amplitude
s = 20*log10( abs( s ) );

% shift time 0 to center of time axis
s = fftshift( s, 3 );

% uncomment this to reuse data
%save wave_fronts_encoded_field_large_scale.mat
%save wave_fronts_encoded_field_small_scale.mat

figure;
imagesc( X, Y, s( :, :, t ) , [ -40 0 ] );
turn_imagesc;
colormap gray;
revert_colormap;
colorbar;
axis square;

xlabel( 'x (m)' );
ylabel( 'y (m)' ); 

graph_defaults;
