function [G, H] = gradhesseval(c,n,w)

% Set up:
  if nargin < 3
      w = 0;
  end
  N = length(c)/2; 
  u = c(1:N); 
  v = c(N+1:end); 
  if ( mod(N, 2) == 0 )
    kk = (-N/2:N/2-1)';
  else
    kk = (-(N-1)/2:(N-1)/2)';
  end
  [t, s] = trigpts(N, [0 2*pi]);

% Initialize G and Hessian H:
  G = zeros(2*N, 1);
  if nargout > 1
    H = zeros(2*N, 2*N);
  end

% Loop over the bodies for the AU term:
  for j = 1:n-1
    a = bsxfun(@times,1-cos(kk'*j*2*pi/n),cos(t*kk')) + ...
      bsxfun(@times,sin(kk'*j*2*pi/n),sin(t*kk'));
    b = bsxfun(@times,-1+cos(kk'*j*2*pi/n),sin(t*kk')) + ...
      bsxfun(@times,sin(kk'*j*2*pi/n),cos(t*kk'));
    f = (a*u + b*v).^2 + (-b*u + a*v).^2;
    dfdu = 2*(bsxfun(@times,a*u + b*v,a) + bsxfun(@times,-b*u + a*v,-b));
    G(1:N) = G(1:N) - n/4*bsxfun(@rdivide,dfdu,f.^(3/2))'*s';
    dfdv = 2*(bsxfun(@times,a*u + b*v,b) + bsxfun(@times,-b*u + a*v,a));
    G(N+1:2*N) = G(N+1:2*N) - n/4*bsxfun(@rdivide,dfdv,f.^(3/2))'*s';
    if nargout > 1
      f52 = f.^(5/2);
      for l = 1:N      
        dfdul = dfdu(:,l);
        dfdvl = dfdv(:,l);
        for k = 1:N
          dfdvk = dfdv(:,k); 
          d2fduldvk = 2*(a(:,l).*b(:,k) - b(:,l).*a(:,k));
          H(l,k+N) = H(l,k+N) - n/4*s*((f.*d2fduldvk-3/2*dfdvk.*dfdul)./f52); 
          if ( k <= l )
            dfduk = dfdu(:,k);
            d2fdulduk = 2*(a(:,l).*a(:,k) + b(:,l).*b(:,k)); 
            H(l,k) = H(l,k) - n/4*s*((f.*d2fdulduk-3/2*dfduk.*dfdul)./f52);
            H(k,l) = H(l,k);
            H(l+N,k+N) = H(l+N,k+N) - ...
              n/4*s*((f.*d2fdulduk-3/2*dfdvk.*dfdvl)./f52);
            H(k+N,l+N) = H(l+N,k+N);
          end  
        end     
      end
    end
  end

% Add the AK term:
  if ( mod(N, 2) == 0 )
    main = 2*pi*n*[w;(kk(2:end)+w);w;(kk(2:end)+w)].^2;
    G = G + main.*[u;v];
    if nargout > 1
      H = H + diag(main);
    end
    idx1 = N/2;
    idx2 = N+N/2;
  else
    main = 2*pi*n*[(kk+w);(kk+w)].^2;
    G = G + main.*[u;v];
    if nargout > 1
      H = H + diag(main);
    end
    idx1 = (N-1)/2;
    idx2 = N+(N-1)/2;
  end

% Get gid of the constant term and use symmetry:
  G = [G(1:idx1); G(idx1+2:idx2); G(idx2+2:end)];
  if nargout > 1
    H = [H(1:idx1,1:idx1), H(1:idx1,idx1+2:idx2), H(1:idx1,idx2+2:end); ...
         H(idx1+2:idx2,1:idx1), H(idx1+2:idx2,idx1+2:idx2), ...
         H(idx1+2:idx2,idx2+2:end); H(idx2+2:end,1:idx1), ...
         H(idx2+2:end,idx1+2:idx2), H(idx2+2:end,idx2+2:end)];
    H(N:end,1:N-1) = H(1:(N-1),N:end)';
  end
  
end
