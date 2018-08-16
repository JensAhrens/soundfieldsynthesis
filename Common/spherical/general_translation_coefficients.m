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

function [ coefficients ] = general_translation_coefficients( ...
                                    N, N_prime, kr, d_alpha, d_beta, type )
%
% This function implements Eq. (E.5) and (E.6), respectively, and
% calculates the general translation coefficients up to n = N-1 and 
% n' = N'-1 .
%
% Indexing of coefficients: 1st dimension ( n , m  ): 1 + n^2  + n  + m
%                           2nd dimension ( n', m' ): 1 + n'^2 + n' + m'
%
% Author: Jens Ahrens

fprintf( 'Calculating general translation coefficients ... ' );

% initialize
coefficients = zeros( N^2, N_prime^2 );

% loop over all possible n, m, n', and m'
for n = 0 : N - 1
    for m = -n : n
        for n_prime = 0 : N_prime - 1
            for m_prime = -n_prime : n_prime

                coefficients( 1 + n^2 + n + m, 1 + n_prime^2 + n_prime + m_prime ) = ...
                            get_coefficient( n, m, n_prime, m_prime, kr, d_alpha, d_beta, type );
            
            end
        end
    end
end

fprintf( 'Done.\n' );

end

function [ coefficient ] = get_coefficient( n, m, n_prime, m_prime, kr, d_alpha, d_beta, type )
% This function actually implements Eq. (E.5) and (E.6), respectively.

if ( strcmp( type, 'EI' ) )
    bessel_hankel = @( n, x ) sphbesselh( n, 2, x );
    
elseif ( strcmp( type, 'II' ) || strcmp( type, 'EE' ) )
    bessel_hankel = @( n, x ) sphbesselj( n, x ); 
    
else
    error( 'Unknown translation type.' );

end

coefficient = 0;

m_prime_prime = m_prime - m;

for n_prime_prime = abs( n_prime - n ) : 2 : n_prime + n

    % only evaluate for non-zero spherical harmonics
    if ( n_prime_prime >= abs( m_prime_prime ) )
       
        coefficient = coefficient + sqrt( ( ( 2*n_prime+1 ) * ( 2*n_prime_prime+1 ) * ( 2*n+1 ) ) / (4*pi) ) .* ...
                1i^( n + n_prime_prime - n_prime ) .* ... (-1)^( n_prime_prime ) .* ...
                    e_symbol( [ m_prime -m -m_prime_prime ], [ n_prime n n_prime_prime ] ) .* ...
                        bessel_hankel( n_prime_prime, kr ) .* ...
                            sphharm( n_prime_prime, m_prime_prime, d_beta, d_alpha );

    end
    
end

end
