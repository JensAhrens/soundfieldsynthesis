function [ sig ] = noise( varargin )

dim1 = varargin{1};
dim2 = varargin{1};
dim3 = 1;

if ( nargin == 1 )
    ;
elseif ( nargin == 2 )
    dim2 = varargin{2};
elseif ( nargin == 3 )
    dim3 = varargin{3};
else 
    error( 'Can only do up to 3 dimensions.' )
end

sig = squeeze( randn( dim1, dim2, dim3 ) );

end

