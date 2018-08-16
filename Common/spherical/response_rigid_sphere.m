function [ S ] = response_rigid_sphere( A, x_mesh, y_mesh, sphere_location, source_type, source_location, k, N )
% source_type: 'plane' or 'point'
%               if 'plane': source_location = [ theta_pw phi_pw ]
%               if 'point': source_location = [ x_s y_s ]
%
% horizontal plane

% translate the coordinate system to consider location of sphere
x_mesh = x_mesh - sphere_location( 1 );
y_mesh = y_mesh - sphere_location( 2 );

% translate the coordinate system to consider location of sphere (it
% doesn't matter for the plane wave)
if ( strcmp( source_type, 'point' ) )
    source_location( 1 ) = source_location( 1 ) - sphere_location( 1 );
    source_location( 2 ) = source_location( 2 ) - sphere_location( 2 );
end    

% coordinate system with origin at sphere
r     = sqrt( x_mesh.^2 + y_mesh.^2 );
alpha = atan2( y_mesh, x_mesh );

% calculate source field
if ( strcmp( source_type, 'plane' ) )
    
    k_x = k .* cos( source_location( 1 ) );
    k_y = k .* sin( source_location( 1 ) );
    
    if ( source_location( 2 ) ~= pi/2 )
        error( 'Horizontal plane only.' );
    end
    
    S = exp( -1i .* ( k_x .* x_mesh + k_y .* y_mesh ) );
    
elseif ( strcmp( source_type, 'point' ) )
    
    r_0 = sqrt( ( x_mesh - source_location( 1 ) ).^2 + ( y_mesh - source_location( 2 ) ).^2 );
    
    S = 1 / (4*pi) .* exp( -1i .* k .* r_0 ) ./ r_0;
    
    % source coordinates relative to sphere
    r_s     = sqrt( source_location( 1 ).^2 + source_location( 2 ).^2 );
    alpha_s = atan2( source_location( 2 ), source_location( 1 ) );
    
end

% for testing
%S = zeros( size( S ) );

S_scat = zeros( size( r ) );

for n = 0 : N 
    %disp( [ 'Calculating order ' num2str( n ) ' of ' num2str( N ) '.' ] ); 
    
    %bessel = sphbesselj( n,    k .* r );
    hankel = sphbesselh( n, 2, k .* r );
    
    bessel_prime = 1/(2*n+1) .* ( n .* sphbesselj( n-1,    k*A ) - (n+1) .* sphbesselj( n+1,    k*A ) ) .* k;
    hankel_prime = 1/(2*n+1) .* ( n .* sphbesselh( n-1, 2, k*A ) - (n+1) .* sphbesselh( n+1, 2, k*A ) ) .* k;
    
    for m = -n : n
        
        % Eq. (2.38)
        if ( strcmp( source_type, 'plane' ) )
            S_breve = 4 .* pi .* 1i^(-n) .* sphharm( n, -m, source_location( 2 ), source_location( 1 ) );
            
        elseif ( strcmp( source_type, 'point' ) )
            S_breve = (-1i) .* k .* sphbesselh( n, 2, k .* r_s ) .* sphharm( n, -m, pi/2, alpha_s );
        else
            error( 'Unknown source_type.' );
        end

        % Eq. (3.108)
        S_breve_scattered = - bessel_prime ./ hankel_prime .* S_breve;

        S_scat = S_scat + S_breve_scattered .* hankel .* sphharm( n, m, pi/2, alpha ); 
        
        % this combines Eq. (2.32a) and (2.32b) 
        % ( incoming sound field plus outgoing scattered sound field )
        %S_int = S_int + ( S_breve .* bessel + S_breve_scattered .* hankel ) .* sphharm( n, m, pi/2, alpha ); 

    end
    
end

S = S + S_scat;

% no sound field inside the sphere
S( r < A ) = 0;

% add loudspeaker signals
%S = sum( S, 2 );

end