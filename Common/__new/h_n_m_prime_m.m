function [ out ] = h_n_m_prime_m( n, m_prime, m, beta )

%%% THIS DOES NOT WORK !!! %%%
% if ( m == 0 )
% 
%     Pnm = legendre( n, cos( beta ) );
% 
%     if n~=0
%         Pnm = squeeze( Pnm ( abs( m_prime ) + 1, :, : ) );
%     end
%     
%     out = (-1)^m_prime * sqrt( factorial( n - abs( m_prime ) ) / factorial( n + abs( m_prime ) ) ) .* Pnm;
%    
%     return;
%     
% end
    
out = 0;

for sigma = max( 0, -(m+m_prime) ) : min( n-m, n-m_prime )

    out = out + ( (-1)^(n-sigma) * cos( beta/2 )^( 2*sigma + m_prime + m ) * sin( beta/2 )^( 2*n - 2*sigma - m_prime - m ) ) / ...
                    ( factorial( sigma ) * factorial( n - m - sigma ) * factorial( n - m_prime - sigma ) * factorial( m + m_prime + sigma ) );

end

out = out * epsilon( m ) * epsilon( m_prime ) * sqrt( factorial( n + m_prime ) * factorial( n - m_prime ) * factorial( n + m ) * factorial( n - m ) );
    
end

function [ out ] = epsilon( m )

if ( m > 0 )
    out = (-1)^m;
else
    out = 1;
end

end