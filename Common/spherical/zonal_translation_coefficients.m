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

function [ coefficients ] = zonal_translation_coefficients( N, N_prime, kr, type )
% Evaluates the zonal translation coefficients from Sec. 3.3.3.1 up to 
% n = N - 1 and n' = N_prime - 1 .
%
% Structure of coefficients: columns: n, rows: n'
%
% type: either 'EI', 'II', or 'EE'
%
% Author: Jens Ahrens

fprintf( 'Calculating zonal translation coefficients ... ' );

% sanity check
if ( ( N == 0 ) || ( N_prime == 0 ) )
    error( 'N and N'' need to be larger than 0.' );
end

% determine which type of translation is asked for
if ( strcmp( type, 'EI' ) )
    bessel_hankel = @( n, x ) sphbesselh( n, 2, x );
    
elseif ( strcmp( type, 'II' ) || strcmp( type, 'EE' ) )
    bessel_hankel = @( n, x ) sphbesselj( n, x ); 
    
else
    error( 'Unknown translation type.' );
    
end

% initialize
coefficients = zeros( 2*N, 2*N_prime );

% calculate initial values for all n with n' = 0, Eq. (3.28)
for n = 0 : 2*N - 1
    coefficients( n+1, 1 ) = (-1)^n * sqrt( 2*n + 1 ) * bessel_hankel( n, kr ); 
end

% calculate initial values for all n' with n = 0, Eq. (3.29)
for n_prime = 0 : 2*N_prime - 1
    coefficients( 1, n_prime+1 ) = sqrt( 2*n_prime + 1 ) * bessel_hankel( n_prime, kr ); 
end

%%% calculate all values of Eq.(3.30) in which the factor a_{-1} arises %%%

for n = 1 : 2*N - 2

    % this is (E|I)_(n,1)
    coefficients( n + 1, 2 ) = 1 / a_n( 0 ) * ( a_n( n - 1 ) * coefficients( n,     1 ) - ...
                                                a_n( n     ) * coefficients( n + 2, 1 ) );
end

for n_prime = 1 : 2*N_prime - 2

    % this is (E|I)_(1,n')
    coefficients( 2, n_prime + 1 ) = 1 / a_n( 0 ) * ( a_n( n_prime - 1 ) * coefficients( 1, n_prime     ) - ...
                                                      a_n( n_prime     ) * coefficients( 1, n_prime + 2 ) );
end

% now fill up the missing values
for n_prime = 1 : 2*N_prime - 1
       
    for n = 2 : 2*N - 2
        % this is (E|I)_(n),(n'+1)
        coefficients( n + 1, n_prime + 2 ) = - 1 / a_n( n_prime ) * ( a_n( n )           * coefficients( n + 2, n_prime + 1 ) - ...
                                                                      a_n( n - 1 )       * coefficients( n    , n_prime + 1 ) - ...
                                                                      a_n( n_prime - 1 ) * coefficients( n + 1, n_prime     ) );
                     
    end
    
end

% remove data that is not required anymore
coefficients = coefficients( 1 : N, 1 : N_prime );

fprintf( 'Done.\n' );

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ a ] = a_n( n )
% this function implements Eq. (3.31)

if ( n < -1 )
    error( 'n has to be larger than -1.' );
elseif ( n == -1 )
    a = 0;
else
    a = ( n+1 ) / sqrt( ( 2*n + 1 ) * ( 2*n + 3 ) );
end

end