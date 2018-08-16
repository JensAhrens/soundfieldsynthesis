%=============================================
%
%   Compute the transfer function
%
%=============================================

function [Spectr hh err] = ComputeTransferFunction(Signal, Reference, frameSize, method)

    if (nargin < 4)
        method = 'spec';
    end

    Signal = Signal(:);
    Reference = Reference(:);

    switch (method)
        case 'spec'
            [Sp h er] = ComputeTransferFunctionSpec(Signal, Reference, frameSize);
            
        case 'corr'
            [Sp h er] = ComputeTransferFunctionCorr(Signal, Reference, frameSize);
            
        case 'optim'
            [Sp h er] = ComputeTransferFunctionOptim(Signal, Reference, frameSize);
            
        case 'mls'
            [Sp h er] = ComputeTransferFunctionMLS(Signal, frameSize);
    end
            
    Spectr = Sp;

    % return impulse response if asked
    if (nargout > 1)
        hh = h;
    end
    
    % return impulse response if asked
    if (nargout > 2)
        err = er;
    end
end


%===================================================
%
%   Compute the transfer function - spectral method
%   tardeoff precision and speed, default
%
%===================================================

function [Spectr hh err] = ComputeTransferFunctionSpec(Signal, Reference, frameSize)

    % convert to frequency domain
    nSamples = int32(length(Signal));

    frameWhole = 2*frameSize;
    ha = sqrt(hann(frameWhole));
    
    begSamples = 1:floor(frameSize/4):nSamples-frameWhole;
    endSamples = begSamples + frameWhole - 1;
    totFrames = length(begSamples);
    
    SigSpec(1:totFrames,1:frameWhole) = 0.0;
    RefSpec(1:totFrames,1:frameWhole) = 0.0;

    for nFrame = 1:totFrames
        SigSpec(nFrame,:) = fft(ha .* Signal(begSamples(nFrame):endSamples(nFrame)));
        RefSpec(nFrame,:) = fft(ha .* Reference(begSamples(nFrame):endSamples(nFrame)));
    end

    % compute the transfer coefficients per bin
    H(1:frameSize) = 0.0;
    for iBin = 1:frameSize
        Input = squeeze(RefSpec(:,iBin));
        Output = squeeze(SigSpec(:,iBin));
        H(iBin) = Input \ Output;
    end
    
    Spectr = H;

    % compute the impulse response
    H(frameSize+1) = 0.0 + i * 0.0;
    H(frameSize+2:2*frameSize) = conj(H(frameSize:-1:2));
    h0 = ifft(H);
    
    % return impulse response if asked
    if (nargout > 1)
        hh = h0; % h0(1:frameSize);
    end
    
    % return error if asked
    if (nargout > 2)
        Out = filter(hh, 1, Reference);
        err = sqrt(mean((Out-Signal).^2));
    end
    
end


%=========================================================
%
%   Compute the transfer function - correlation method
%   fastest and least precise
%
%=========================================================

function [Spectr hh err] = ComputeTransferFunctionCorr(Signal, Reference, frameSize)

    % estimate the trasfer function
    [CC lags] = xcorr(Signal, Reference);
    [val index0] = min(abs(lags));
    h0 = CC(index0:index0+2*frameSize-1);

    % scale it
    RR = ComputeBestRatio(Signal, Reference, h0);
    h0 = RR * h0;

    % compute the spectrum
    H = fft(h0);
    Spectr = H(1:frameSize).';
    
    % return impulse response if asked
    if (nargout > 1)
        hh = h0(1:frameSize);
    end
    
    % return error if asked
    if (nargout > 2)
        Out = filter(hh, 1, Reference);
        err = sqrt(mean((Out-Signal).^2));
    end
end


%===================================================
%
%   Compute Trasfer Function - optimization method
%   most precise, but much, much slower
%
%===================================================

function [Spectr hh err] = ComputeTransferFunctionOptim(Signal, Reference, frameSize)

    eps = 1e-5;

    %
    %   In-line function for optimization criterion
    %
    function [res gradX] = crit(x)

        Out = filter(x, 1, Reference);
        res = sqrt(mean((Out-Signal).^2));

        % compute the derivatives if asked
        if (nargout > 1)
            curFun = res;
            numberOfVariables = length(x);
            gradX(1:numberOfVariables) = 0.0;
            for iPar = 1:numberOfVariables
                grX = x;
                grX(iPar) = grX(iPar) + eps;
                grFun = crit(grX);
                gradX(iPar) = (grFun - curFun) / eps;
            end
        end
    end

    %
    %   Find the best filter 
    %
    [Spectr h0] = ComputeTransferFunctionSpec(Signal, Reference, frameSize);

    options = fminunc('defaults');
    options = optimset(options,'Display','notify', 'TolX',eps, 'TolFun',1000*eps, 'MaxFunEvals', 100, 'GradObj','on', 'LargeScale', 'off');
    h = fminunc(@crit, h0, options); 

    H = fft(h);
    Spectr = H(1:frameSize).';

    if (nargout > 1)
        hh = h;
    end

    % return error if asked
    if (nargout > 2)
        Out = filter(hh, 1, Reference);
        err = sqrt(mean((Out-Signal).^2));
    end
end


%=============================================
%
%   Compute Best Ratio - helper function
%
%=============================================

function Ratio = ComputeBestRatio(Signal, Reference, hh)

    eps = 1e-5;

    %
    %   In-line function for optimization criterion
    %
    function [res gradX] = crit(x)

        Out = filter(x(2)*hh, 1, Reference);
        res = sqrt(mean((Out-Signal-x(1)).^2));

        % compute the derivatives if asked
        if (nargout > 1)
            curFun = res;
            numberOfVariables = length(x);
            gradX(1:numberOfVariables) = 0.0;
            for iPar = 1:numberOfVariables
                grX = x;
                grX(iPar) = grX(iPar) + eps;
                grFun = crit(grX);
                gradX(iPar) = (grFun - curFun) / eps;
            end
        end
    end

    %
    %   Find the best ratio 
    %
    r0 = [0.0 1.0];
    
    options = fminunc('defaults');
    options = optimset(options,'Display','notify', 'TolX',eps, 'TolFun',eps, 'MaxFunEvals', 100, 'GradObj','on', 'LargeScale', 'off');
    RR = fminunc(@crit, r0, options); 
    
    Ratio = RR(2);

end


%===================================================
%
%   Compute Trasfer Function - MLS Method
%   Precise and fast but requires structured reference
%
%===================================================

function [Spectr hh err] = ComputeTransferFunctionMLS(Signal, frameSize)
    Spectr = [];
    err=[];
    
    if(isstruct(framesize))
        offset=frameSize.offset;
        reps=frameSize.reps;
        N=frameSize.N;
        DCCoupling=frameSize.N;
        useimp=frameSize.useimp;
    else
        offset=0;
        reps=1;
        N=11;
        DCCoupling=0;
        useimp=0;
    end
    hh= AnalyseMLSSequence(Signal,offset,reps,N,DCCoupling,useimp);
end


function [ out ] = H_n_m_prime_m( n, m_prime, m, beta )

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

out = out * Epsilon( m ) * Epsilon( m_prime ) * sqrt( factorial( n + m_prime ) * factorial( n - m_prime ) * factorial( n + m ) * factorial( n - m ) );
    
end

function [ out ] = Epsilon( m )

if ( m > 0 )
    out = (-1)^m;
else
    out = 1;
end

end
