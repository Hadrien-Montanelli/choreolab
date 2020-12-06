function [A,G] = actiongradeval(c,n,w)

% Set up:
  if nargin < 3
      w = 0;
  end
  N = length(c)/2; 
  u = c(1:N); 
  v = c(N+1:end);
  c = u + 1i*v;
  if ( mod(N, 2) == 0 )
    k = (-N/2:N/2-1)';
  else
    k = (-(N-1)/2:(N-1)/2)';
  end
  [t,s] = trigpts(N,[0 2*pi]);

% Evaluate action A:
  c = bsxfun(@times,exp(1i*k*(0:n-1)*2*pi/n),c(:,1));
  vals = ifft(N*c);
  dc = 1i*k.*c(:,1);
  if ( mod(N, 2) == 0 )
    dc(1,1) = 0; 
  end
  dvals = ifft(N*dc);
  K = abs(dvals+1i*w*vals(:,1)).^2; 
  U = 0;
  for i = 2:n
    U = U - 1./abs(vals(:,1)-vals(:,i)); 
  end
  A = n/2*s*(K-U);
  
if nargout > 1
    
% Initialize gradient G:
  G = zeros(2*N,1);
  
% Loop over the bodies for the AU contribution:
  for j = 1:n-1
    a = bsxfun(@times,1-cos(k'*j*2*pi/n),cos(t*k')) + ...
      bsxfun(@times,sin(k'*j*2*pi/n),sin(t*k'));
    b = bsxfun(@times,-1+cos(k'*j*2*pi/n),sin(t*k')) + ...
      bsxfun(@times,sin(k'*j*2*pi/n),cos(t*k'));
    f = (a*u + b*v).^2 + (-b*u + a*v).^2;
    df = 2*(bsxfun(@times,a*u + b*v,a) + bsxfun(@times,-b*u + a*v,-b));
    G(1:N) = G(1:N) - n/4*bsxfun(@rdivide,df,f.^(3/2))'*s';
    df = 2*(bsxfun(@times,a*u + b*v, b) + bsxfun(@times,-b*u + a*v,a));
    G(N+1:2*N) = G(N+1:2*N) - n/4*bsxfun(@rdivide,df,f.^(3/2))'*s';
  end
  
% Add the AK contribution:
  if ( mod(N, 2) == 0 )
    G = G + 2*pi*n*[w;(k(2:end)+w);w;(k(2:end)+w)].^2.*[u;v];
  else
    G = G + 2*pi*n*[(k+w);(k+w)].^2.*[u;v];
  end
end

end