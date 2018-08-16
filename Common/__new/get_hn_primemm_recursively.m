function [ H ] = get_hn_primemm_recursively( N, beta, factorials )
% N          : order (n is evaluated from 0 : N)
% beta       : Euler rotation angle beta_E in rad
% factorials : [ 1 x 4*N-4 ] vector with precalculated values of the 
%              factorial from 0 to 4*N-3
%
% H          : rotation coefficients H_n^(m',m)
%              indexing: [ n, m', m ]
%              size: [ N+1, 2N+1, 2N+1 ], not all entries carry information

beta = mod( beta + pi, 2*pi ) - pi;

P = 2*N;

if ( beta < 0 )
    betaIsNegative = true;
    beta = abs( beta );
else
    betaIsNegative = false;
end

H = zeros( P, 2*P+1, 2*P+1 );

% initial values
for n = 0 : P 
    for m_prime = -n : n
        H( n+1, m_prime+P+1, P+1 ) = (-1)^m_prime * sqrt( factorials( n - abs( m_prime ) + 1 ) / factorials( n + abs( m_prime ) + 1 ) ) * AssLegendre( n, abs( m_prime ), cos( beta ) );         
    end
end

% recursive relation
for m = 0 : P-2
    for n = max( 2, m ) : P
        for m_prime = -n+1 : n-1
            
            if ( n == m || n == m+1 )
                H( n, m_prime+P+1, m+P+2 ) = 0;
                
            else
                H( n, m_prime+P+1, m+P+2 ) = 1 / Bnm( n, m ) * ...
                    ( .5 * ( Bnm( n, -m_prime-1 ) * ( 1-cos( beta ) ) * H( n+1, m_prime+P+2, m+P+1 ) - ...
                             Bnm( n,  m_prime-1 ) * ( 1+cos( beta ) ) * H( n+1, m_prime+P, m+P+1 ) ) - ... 
                             Anm( n-1, m_prime ) * sin( beta ) * H( n+1, m_prime+P+1, m+P+1 ) );
            end

        end
    end
end

% remove faulty overhead
H = H( 1 : N+1, (end+1)/2-N : (end+1)/2+N, (end+1)/2-N : (end+1)/2+N );

if ( betaIsNegative )
    for n = 0 : N
        for m = 0 : n
            for m_prime = -n : n
                H( n+1, m_prime+N+1, m+N+1 ) = (-1)^( m + m_prime ) * H( n+1, m_prime+N+1, m+N+1 );
            end
        end
    end
end

% symmetries
for n = 0 : N
	for m = -n : -1
        for m_prime = 0 : n
            H( n+1, m_prime+N+1, m+N+1 ) = H( n+1, m+N+1, m_prime+N+1 );
        end
        for m_prime = -n : -1
            H( n+1, m_prime+N+1, m+N+1 ) = H( n+1, -m_prime+N+1, -m+N+1 );
        end
	end
end

end

function [ a ] = Anm( n, m )

m = abs( m );

if ( m <= n )
    a = sqrt( ( (n+1+m)*(n+1-m) ) / ( (2*n+1)*(2*n+3) ) );
else 
    a = 0;
end

end

function [ b ] = Bnm( n, m )

if ( 0 <= m )
    b = sqrt( ( (n-m-1)*(n-m) ) / ( (2*n-1)*(2*n+1) ) );
elseif ( -n <= m && m < 0 )
    b = - sqrt( ( (n-m-1)*(n-m) ) / ( (2*n-1)*(2*n+1) ) );
elseif ( abs( m ) > n )
    b = 0;
    warning( 'abs( m ) > n. This should not happen.' );
else
    b = 0;
    warning( 'This line may never be reached.' );
end

end
