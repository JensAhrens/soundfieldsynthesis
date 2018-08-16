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

function [ coefficients ] = sectorial_translation_coefficients( M, N_prime, kr, d_alpha, d_beta, type )
% Evaluates the sectorial translation coefficients from Sec. 3.5.3.1 up to 
% m = M and n' = N_prime - 1 .
%
% Indexing of coefficients: 1st dimension ( n', m' ): 1 + n'^2 + n' + m'
%                           2nd dimension:  m
%
% type: either 'EI', 'II', or 'EE'
%
% Author: Jens Ahrens

fprintf( 'Calculating sectorial translation coefficients ... ' );

% sanity check
if ( ( M == 0 ) || ( N_prime == 0 ) )
    error( 'M and N'' need to be larger than 0.' );
end

% Eq. (3.55)
coefficients_3_55 = sectorial_translation_coefficients_2( N_prime, M, kr, d_alpha, d_beta, type );

% initialize
n_primes = zeros( 1, 0 );

for n_prime = 0 : N_prime - 1
    n_primes = [ n_primes, repmat( n_prime, [ 1, 2*n_prime + 1 ] ) ];
end

[ ms, n_primes ] = meshgrid( -M : M, n_primes );

% from Eq. (3.55)
coefficients = (-1).^( abs( ms ) + n_primes ) .* coefficients_3_55;

fprintf( 'Done.\n' );

end

function [ coefficients ] = sectorial_translation_coefficients_2( N_prime, M, kr, d_alpha, d_beta, type )
% evaluates the coefficients on the right hand side of Eq. (3.55)

% determine which type of translation is asked for
if ( strcmp( type, 'EI' ) )
    bessel_hankel = @( n, x ) sphbesselh( n, 2, x );
    
elseif ( strcmp( type, 'II' ) || strcmp( type, 'EE' ) )
    bessel_hankel = @( n, x ) sphbesselj( n, x ); 
    
else
    error( 'Unknown translation type.' );
    
end

% 1st dimension:  1 + n'^2 + n' + m'
% 2nd dimension:  m
% Put enough headroom.
coefficients = zeros( ( 2*N_prime )^2, 4*M + 1 );

m = 0;

% Eq. (3.56)
for n_prime = 0 : 2*N_prime - 1
    for m_prime = -n_prime : n_prime

        coefficients( 1 + n_prime^2 + n_prime + m_prime, m + 2*M + 1 ) = ...
                                sqrt( 4*pi ) * (-1)^n_prime * bessel_hankel( n_prime, kr ) * ...
                                    sphharm( n_prime, -m_prime, d_beta, d_alpha ); 
   
    end
end

n_prime = 0;
m_prime = 0;

% Eq. (3.57)
for m = -2*M : 2*M
    coefficients( 1 + n_prime^2 + n_prime + m_prime, m + 2*M + 1 ) = ...
                            sqrt( 4*pi ) * bessel_hankel( abs( m ), kr ) * ...
                                sphharm( abs( m ), m, d_beta, d_alpha );     
    
end

% This loop evaluates Eq. (E.7) and (E.8).
for m = 0 : 2*M - 1
    for n_prime = 1 : 2*N_prime - 2
        for m_prime = -n_prime : n_prime
            
            % Eq. (E.7)
            coefficients( 1 + n_prime^2 + n_prime + m_prime, - m + 2*M ) = 1 / b_nm( m+1, -m-1 ) .* ...
                                ( b_nm( n_prime    ,   m_prime     ) * coefficients( 1 + ( n_prime-1 )^2 + ( n_prime-1 ) + m_prime + 1, - m + 2*M + 1 ) - ...
                                  b_nm( n_prime + 1, - m_prime - 1 ) * coefficients( 1 + ( n_prime+1 )^2 + ( n_prime+1 ) + m_prime + 1, - m + 2*M + 1 ) ); 

            % Eq. (E.8)
            if ( ( 1 + ( n_prime-1 )^2 + ( n_prime-1 ) + m_prime - 1 ) < 1 )       
                
                coefficients( 1 + n_prime^2 + n_prime + m_prime, m + 2*M + 2 ) = 1 / b_nm( m+1, -m-1 ) .* ...
                                ( - b_nm( n_prime + 1,   m_prime - 1 ) * coefficients( 1 + ( n_prime+1 )^2 + ( n_prime+1 ) + m_prime - 1, m + 2*M + 1 ) );
            
            else
            
                coefficients( 1 + n_prime^2 + n_prime + m_prime, m + 2*M + 2 ) = 1 / b_nm( m+1, -m-1 ) .* ...
                                ( b_nm( n_prime    , - m_prime     ) * coefficients( 1 + ( n_prime-1 )^2 + ( n_prime-1 ) + m_prime - 1, m + 2*M + 1 ) - ...
                                  b_nm( n_prime + 1,   m_prime - 1 ) * coefficients( 1 + ( n_prime+1 )^2 + ( n_prime+1 ) + m_prime - 1, m + 2*M + 1 ) );

            end
            
        end     
    end
end

% Finally, remove unrequired data.
coefficients = coefficients( 1 : N_prime^2, M + 1 : end - M );

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ b ] = b_nm( n, m )
% This function implements Eq. (3.31).

if ( ( 0 <= m ) && ( m <= n ) )
    sign = 1;

elseif ( ( -n <= m ) && ( m <= 0 ) )
    sign = -1;

elseif ( abs( m ) > n )
    b = 0;
    return;

else 
    error( 'b\_nm is undefined for this combination of n and m.' );

end

b = sign * sqrt( ( ( n-m-1 ) * ( n-m ) ) / ( ( 2*n-1 ) * ( 2*n+1 ) ) );

end
