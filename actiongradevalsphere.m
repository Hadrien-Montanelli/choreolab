function [A, G] = actiongradevalsphere(c,n,w,R)

% Set up:
  if nargin < 3
    w = 0;
    R = 2;
  end
  if nargin < 4
    R = 2;
  end
  N = length(c)/2; 
  u = c(1:N); 
  v = c(N+1:end);
  c = u + 1i*v;
  if ( mod(N, 2) == 0 )
    k = (-N/2:N/2-1)';
    kz = [0,-N/2+1:N/2-1]';
  else
    k = (-(N-1)/2:(N-1)/2)';
    kz = k;
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
  v1 = abs(vals(:,1));
  U = 0;
  for j = 2:n
    vj = abs(vals(:,j));
    D = 2*R^2*abs(vals(:,1)-vals(:,j))./(sqrt(R^2+v1.^2).*sqrt(R^2+vj.^2));
    U = U - 1/R*(2*R^2-D.^2)./(D.*sqrt(4*R^2-D.^2));
  end
  K = (2*R^2*abs(dvals+1i*w*vals(:,1))./(R^2+v1.^2)).^2;
  A = s*n/2*(K-U);
  
if nargout > 1
% Initialize gradient G:
  G = zeros(2*N,1);

% Loop over the bodies for the AU contribution:
  a0 = cos(t*k');
  b0 = -sin(t*k');
  r0 = sqrt(R^2 + (a0*u + b0*v).^2 + (-b0*u + a0*v).^2);
  dr0du = bsxfun(@times,a0*u + b0*v,a0) + bsxfun(@times,-b0*u + a0*v,-b0);
  dr0du = bsxfun(@rdivide,dr0du,r0);
  dr0dv = bsxfun(@times,a0*u + b0*v,b0) + bsxfun(@times,-b0*u + a0*v,a0);
  dr0dv = bsxfun(@rdivide,dr0dv,r0);
  for j = 1:n-1
    c = bsxfun(@times,cos(k'*j*2*pi/n),cos(t*k')) - ...
      bsxfun(@times,sin(k'*j*2*pi/n),sin(t*k'));
    d = -bsxfun(@times,cos(k'*j*2*pi/n),sin(t*k')) - ...
      bsxfun(@times,sin(k'*j*2*pi/n),cos(t*k'));
    r = sqrt(R^2 + (c*u + d*v).^2 + (-d*u + d*v).^2);
    a = a0 - c;
    b = b0 - d;
    f = (a*u + b*v).^2 + (-b*u + a*v).^2;
    D = 2*R^2*sqrt(f)./(r0.*r);
    dr = bsxfun(@times,c*u + d*v,c) + bsxfun(@times,-d*u + c*v,-d);
    dr = bsxfun(@rdivide,dr,r);
    df = 2*(bsxfun(@times,a*u + b*v,a) + bsxfun(@times,-b*u + a*v,-b));
    dD = 2*R^2*(bsxfun(@rdivide,df,r0.*r.*(2*sqrt(f))) - ...
      bsxfun(@rdivide,dr0du,r0.^2.*r./sqrt(f)) - ...
      bsxfun(@rdivide,dr,r0.*r.^2./sqrt(f)));
    G(1:N) = G(1:N) + ...
      n/(2*R)*bsxfun(@rdivide,-8*R^4*dD,D.^2.*(4*R^2-D.^2).^(3/2))'*s';
    dr = bsxfun(@times,c*u + d*v,d) + bsxfun(@times,-d*u + c*v,c);
    dr = bsxfun(@rdivide,dr,r);
    df = 2*(bsxfun(@times,a*u + b*v,b) + bsxfun(@times,-b*u + a*v,a));
    dD = 2*R^2*(bsxfun(@rdivide,df,r0.*r.*(2*sqrt(f))) - ...
      bsxfun(@rdivide,dr0dv,r0.^2.*r./sqrt(f)) - ...
      bsxfun(@rdivide,dr,r0.*r.^2./sqrt(f)));
    G(N+1:2*N) = G(N+1:2*N) + ...
      n/(2*R)*bsxfun(@rdivide,-8*R^4*dD,D.^2.*(4*R^2-D.^2).^(3/2))'*s';
  end
  
% Add the AK contribution:
  g1 = bsxfun(@times,kz'+w,cos(t*k'));
  g2 = bsxfun(@times,kz'+w,sin(t*k'));
  g = (g1*u - g2*v).^2 + (g2*u + g1*v).^2;
  h = (R^2 + (a0*u + b0*v).^2 + (-b0*u + a0*v).^2).^2;
  dg = 2*(bsxfun(@times,g1*u - g2*v,g1) + bsxfun(@times,g2*u + g1*v,g2));
  dh = (bsxfun(@times,a0*u + b0*v,a0) + bsxfun(@times,-b0*u + a0*v,-b0));
  dh = bsxfun(@times,dh,4*sqrt(h));
  G(1:N) = G(1:N) + 2*R^4*n*(bsxfun(@rdivide,dg,h)'*s' - ...
    bsxfun(@rdivide,dh,h.^2./g)'*s');
  dg = 2*(bsxfun(@times,g2*u + g1*v,g1) + bsxfun(@times,g1*u - g2*v,-g2));
  dh = (bsxfun(@times,a0*u + b0*v,b0) + bsxfun(@times,-b0*u + a0*v,a0));
  dh = bsxfun(@times,dh,4*sqrt(h));
  G(N+1:2*N) = G(N+1:2*N) + 2*R^4*n*(bsxfun(@rdivide,dg,h)'*s' - ...
    bsxfun(@rdivide,dh,h.^2./g)'*s');
end

end