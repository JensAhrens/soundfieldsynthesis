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

% Note that the figures created by this script will never look exactly like
% Fig. 5.24 and 5.25 because of the random components that are involved.

figure_number = '5.24'; % for Fig. 5.24(a) and (b)
%figure_number = '5.25'; % for Fig. 5.25(a) and (b)

if ( strcmp( figure_number, '5.24' ) )
    L = 4;
elseif ( strcmp( figure_number, '5.25' ) )
    L = 0.5;
end

fs = 44100; % temporal sampling frequency
f  = linspace( 0, fs/2, 601 ).';

omega = 2.*pi.*f;
c     = 343;
k     = omega./c;

no_of_simultaneous_vibration_modes = 2;

ratio = 2.5; % section length / lambda, refer to p. 205 bottom
y     = 3; 

frequency_bands_Hz   = [ 0 180 300 420 540 660 780 900, ...
                         1075 1225 1375 1525 1700 1900 2100 2300 2550 2850 3150, ... 
                         3475 4000 5200 6760 8788 11424 14852 ];

frequency_bands_bins = floor( frequency_bands_Hz / (fs/2) * length( f ) ) + 1;

frequency_bands_Hz   = [ frequency_bands_Hz  , f( end ) ];
frequency_bands_bins = [ frequency_bands_bins, length( f ) ];


% initialize etas and weights
etas    = zeros( length( f ), no_of_simultaneous_vibration_modes );
weights = ones(  length( f ), no_of_simultaneous_vibration_modes ) ./ no_of_simultaneous_vibration_modes;
        
% loop over frequency_bands and create random vibration modes
for index_bands = 1 : length( frequency_bands_bins ) - 1
    for index_eta = 1 : no_of_simultaneous_vibration_modes
        
        f_prime = frequency_bands_Hz( index_bands + 1 );
        
        eta_prime = round( ( L*f_prime ) / ( ratio * c ) - 1 );
        
        if ( index_eta ~= 1 )
           eta_prime = abs( round( 2 * randn(1) * eta_prime ) ) + 1;
        end
        
        % save current eta
        etas( frequency_bands_bins( index_bands ) : ...
            frequency_bands_bins(index_bands+1) - 1, index_eta ) = eta_prime;
        
        % save current weight    
        weight_prime = ( randn( 1 ) + i * randn( 1 ) );        
            
        weights( frequency_bands_bins( index_bands ) : ...
        frequency_bands_bins( index_bands + 1 ) - 1, index_eta ) = weight_prime; 
    
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%% prepare spatial fft %%%%%%%%%%%%%%%%%%%%%%%%%%%%
spatial_interval = [ -10 10 ];

delta_x = .001; % sampling interval for spatial fft in meters

X = spatial_interval( 1 ) : delta_x : spatial_interval( 2 );

k_x_s = (2*pi) / delta_x; % spatial sampling frequency

% create k_x
k_x    = linspace( 0, k_x_s/2, ( length( X ) + 1 ) / 2 ); % positive frequencies
k_x(1) = k_x(2); % to avoid numerical instabilities
k_x    = [ -fliplr( k_x( 2 : end ) ), k_x ]; % adds negative frequencies

[ k_x_m, f_m ] = meshgrid( k_x, f );

clear f_m; % we don't need it...
%%%%%%%%%%%%%%%%%%%%%%% end prepare spatial fft %%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize w_tilde
w_tilde = zeros( size( k_x_m ) );
     
% this loop evaluates the components of Eq. (5.24) over all frequencies
for bin = 1 : size( f, 1 )

    disp( [ 'Calculating frequency bin ' num2str( bin ) ' of ' num2str( size( f, 1 ) ) '.' ] );

    % loop over different eta for current frequency band
    for index_eta = 1 : no_of_simultaneous_vibration_modes
        
        eta     = etas(    bin, index_eta );
        weight  = weights( bin, index_eta );

        % loop over vibrating sections
        for l = 0 : eta

            % Eq. (5.21)
            x_2 = -L/2 + ( l + 1 ) * L / ( eta + 1 );
            x_1 = -L/2 + ( l     ) * L / ( eta + 1 );

            % Eq. (5.25), first case
            w_tilde( bin, k_x == 0 ) = w_tilde( bin, k_x == 0 ) + ...
                 weight .* (-1)^l .* ( x_2 - x_1 );
             
            % Eq. (5.25), second case
            w_tilde( bin, k_x ~= 0 ) = w_tilde( bin, k_x ~= 0 ) + ...
                 weight .* (-1)^l .* ( exp( 1i .* k_x( k_x ~= 0 ) .* x_2 ) - exp( 1i .* k_x( k_x ~= 0 ) .* x_1 ) );

        end
    end
    
    % factor 1/(i k_x) from Eq. (5.25), second case
    w_tilde( bin, k_x ~= 0 ) = w_tilde( bin, k_x ~= 0 ) ./ ( 1i .* k_x( k_x ~= 0 ) );

    % Eq. (C.10), first case
    G_kx( bin, abs( k_x ) <= k( bin ) ) = -i/4 * ...
        besselh( 0, 2, sqrt( k( bin ).^2 - k_x( abs( k_x ) <= k( bin ) ).^2 ) .* y );

    % Eq. (C.10), second case
    G_kx( bin, abs( k_x ) > k( bin ) ) = 1/(2*pi) * ...
        besselk( 0, sqrt( k_x( abs( k_x ) > k( bin ) ).^2 - k( bin ).^2 ) .* y );

end

% correct numerical instabilities
w_tilde( isnan( w_tilde ) ) = 0;   
G_kx(    isnan( G_kx    ) ) = 0;

% Eq. (5.26)
S_kx = w_tilde .* G_kx;

% clear some memory
clear w_tilde G_kx k_x_m f_m;

S = ifftx( S_kx, [], 2 );

% clear some memory
clear S_kx;

% pick ear positions, avoid symmetry 
x_left_ear  =  0.09;
x_right_ear = -0.07;

transfer_function = [ S( :, (end-1)/2 + x_left_ear  / delta_x ), ...
                      S( :, (end-1)/2 + x_right_ear / delta_x ) ];

% apply factor in front of integral in (5.20) ( skip the normalization factors... )
transfer_function = transfer_function .* [ sqrt( 1i .* f ), sqrt( 1i .* f ) ];

% numerical Fourier transform
d = real( ifft( [ transfer_function; conj( flipud( transfer_function( 2 : end, : ) ) ) ] ) );

% normalize
d = d ./ max( abs( d( : ) ) );

% Eq. (5.34)
[ coherence, freq ] = mscohere( d( :, 1 ), d( :, 2 ) );

figure;
semilogx( freq ./ pi .* fs/2, coherence, 'k', 'LineWidth', 2 );
axis tight;
grid on
xlabel( 'f (Hz)' );
ylim( [ 0 1 ] );
graph_defaults;

figure;
t = ( 1 : size( d, 1 ) ) ./ fs * 1000; % ms
% left ear (the circshift is necessary to compensate for time aliasing)
plot( t, circshift( d( :, 1 ), [ -60 0 ] ), 'Color', [ .7 .7 .7 ], 'LineWidth', 2 );
hold on;
% right ear (the circshift is necessary to compensate for time aliasing)
plot( t, circshift( d( :, 2 ), [ -60 0 ] ), 'k' );
hold off;
grid on;
xlabel( 't (ms)' );
ylim( [ -1  1 ] );
xlim( [  0 20 ] );
graph_defaults;
