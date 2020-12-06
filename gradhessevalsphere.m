function [G, H] = gradhessevalsphere(c,n,w,R)

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
  if ( mod(N, 2) == 0 )
    kk = (-N/2:N/2-1)';
    kz = [0,-N/2+1:N/2-1]';
  else
    kk = (-(N-1)/2:(N-1)/2)';
    kz = kk;
  end
  [t,s] = trigpts(N,[0 2*pi]);

% Initialize gradient and Hessian:
  G = zeros(2*N, 1);
  if nargout > 1
    H = zeros(2*N, 2*N);
  end
 
% Loop over the bodies for the AU contribution:
  a0 = cos(t*kk');
  b0 = -sin(t*kk');
  r0 = sqrt(R^2 + (a0*u + b0*v).^2 + (-b0*u + a0*v).^2);
  dr0du = bsxfun(@times,a0*u + b0*v,a0) + bsxfun(@times,-b0*u + a0*v,-b0);
  dr0du = bsxfun(@rdivide,dr0du,r0);
  dr0dv = bsxfun(@times,a0*u + b0*v,b0) + bsxfun(@times,-b0*u + a0*v,a0);
  dr0dv = bsxfun(@rdivide,dr0dv,r0);
  for j = 1:n-1
    c = bsxfun(@times,cos(kk'*j*2*pi/n),cos(t*kk')) - ...
      bsxfun(@times,sin(kk'*j*2*pi/n),sin(t*kk'));
    d = -bsxfun(@times,cos(kk'*j*2*pi/n),sin(t*kk')) - ...
      bsxfun(@times,sin(kk'*j*2*pi/n),cos(t*kk'));
    r = sqrt(R^2 + (c*u + d*v).^2 + (-d*u + c*v).^2);
    a = a0 - c;
    b = b0 - d;
    f = (a*u + b*v).^2 + (-b*u + a*v).^2;
    D = 2*R^2*sqrt(f)./(r0.*r);
    drdu = bsxfun(@times,c*u + d*v,c) + bsxfun(@times,-d*u + c*v,-d);
    drdu = bsxfun(@rdivide,drdu,r);
    dfdu = 2*(bsxfun(@times,a*u + b*v,a) + bsxfun(@times,-b*u + a*v,-b));
    dDdu = 2*R^2*(bsxfun(@rdivide,dfdu,r0.*r.*(2*sqrt(f))) - ...
      bsxfun(@rdivide,dr0du,r0.^2.*r./sqrt(f)) - ...
      bsxfun(@rdivide,drdu,r0.*r.^2./sqrt(f)));
    G(1:N) = G(1:N) + ...
        n/(2*R)*bsxfun(@rdivide,-8*R^4*dDdu,D.^2.*(4*R^2-D.^2).^(3/2))'*s';
    drdv = bsxfun(@times,c*u + d*v,d) + bsxfun(@times,-d*u + c*v,c);
    drdv = bsxfun(@rdivide,drdv,r);
    dfdv = 2*(bsxfun(@times,a*u + b*v,b) + bsxfun(@times,-b*u + a*v,a));
    dDdv = 2*R^2*(bsxfun(@rdivide,dfdv,r0.*r.*(2*sqrt(f))) - ...
      bsxfun(@rdivide,dr0dv,r0.^2.*r./sqrt(f)) - ...
      bsxfun(@rdivide,drdv,r0.*r.^2./sqrt(f)));
    G(N+1:2*N) = G(N+1:2*N) + ...
        n/(2*R)*bsxfun(@rdivide,-8*R^4*dDdv,D.^2.*(4*R^2-D.^2).^(3/2))'*s';
    if nargout > 1
      denom0 = ((r0.*r).^3.*f.^(3/2));
      denom1 = (D.^3.*(4*R^2-D.^2).^(3/2));
      denom2 = (D.*(4*R^2-D.^2).^(5/2));
      denom3 = (D.^2.*(4*R^2-D.^2).^(3/2));
      for l = 1:N
        dr0dul = dr0du(:,l);
        dr0dvl = dr0dv(:,l);
        drdul = drdu(:,l);
        drdvl = drdv(:,l);
        dfdul = dfdu(:,l);
        dfdvl = dfdv(:,l);
        dDdul = dDdu(:,l);
        dDdvl = dDdv(:,l);
        for k = 1:N
          dr0duk = dr0du(:,k);
          dr0dvk = dr0dv(:,k);
          drduk = drdu(:,k);
          drdvk = drdv(:,k);
          dfduk = dfdu(:,k);
          dfdvk = dfdv(:,k);
          dDduk = dDdu(:,k);
          dDdvk = dDdv(:,k);
          q = dfdvk/2.*r0.*r - f.*r.*dr0dvk - f.*r0.*drdvk;
          d2r0duldvk = (b0(:,k).*a0(:,l)-a0(:,k).*b0(:,l)-dr0dul.*dr0dvk)./r0;
          d2rduldvk = (d(:,k).*c(:,l)-c(:,k).*d(:,l)-drdul.*drdvk)./r;
          d2fduldvk = 2*(b(:,k).*a(:,l) - a(:,k).*b(:,l));
          dqdul = .5*(d2fduldvk.*r0.*r + dfdvk.*dr0dul.*r + dfdvk.*r0.*drdul)...
            - dfdul.*r.*dr0dvk - f.*drdul.*dr0dvk - f.*r.*d2r0duldvk ...
            - dfdul.*r0.*drdvk - f.*dr0dul.*drdvk - f.*r0.*d2rduldvk;
          d2Dduldvk = 2*R^2*(dqdul.*r0.*r.*f - ...
            q.*(2*r0.*f.*drdul + 2.*r.*f.*dr0dul + r0.*r.*dfdul/2))./denom0; 
          H(l,k+N) = H(l,k+N) + n/(2*R)*s*(16*R^4*(dDdul.*dDdvk)./denom1 - ...
              24*R^4*(dDdul.*dDdvk)./denom2 - 8*R^4*(d2Dduldvk)./denom3);
          if k <= l  
            d2r0dvldvk = (a0(:,k).*a0(:,l)+b0(:,k).*b0(:,l)-dr0dvk.*dr0dvl)./r0;
            d2rdvldvk = (c(:,k).*c(:,l)+d(:,k).*d(:,l)-drdvk.*drdvl)./r;
            d2fdvldvk = 2*(a(:,k).*a(:,l) + b(:,k).*b(:,l));
            dqdvl = .5*(d2fdvldvk.*r0.*r + dfdvk.*dr0dvl.*r + dfdvk.*r0.*drdvl) ...
              - dfdvl.*r.*dr0dvk - f.*drdvl.*dr0dvk - f.*r.*d2r0dvldvk ...
              - dfdvl.*r0.*drdvk - f.*dr0dvl.*drdvk - f.*r0.*d2rdvldvk;
            d2Ddvldvk = 2*R^2*(dqdvl.*r0.*r.*f - ...
              q.*(2*r0.*f.*drdvl + 2.*r.*f.*dr0dvl + r0.*r.*dfdvl/2))./denom0; 
            H(l+N,k+N) = H(l+N,k+N) + ...
              n/(2*R)*s*(16*R^4*(dDdvl.*dDdvk)./denom1 - ...
              24*R^4*(dDdvl.*dDdvk)./denom2 - 8*R^4*(d2Ddvldvk)./denom3); 
            H(k+N,l+N) = H(l+N,k+N);
            q = dfduk/2.*r0.*r - f.*r.*dr0duk - f.*r0.*drduk;
            d2r0dulduk = (a0(:,k).*a0(:,l)+b0(:,k).*b0(:,l)-dr0duk.*dr0dul)./r0;
            d2rdulduk = (c(:,k).*c(:,l)+d(:,k).*d(:,l)-drduk.*drdul)./r;
            dqdul = .5*(d2fdvldvk.*r0.*r + dfduk.*dr0dul.*r + dfduk.*r0.*drdul) ...
              - dfdul.*r.*dr0duk - f.*drdul.*dr0duk - f.*r.*d2r0dulduk ...
              - dfdul.*r0.*drduk - f.*dr0dul.*drduk - f.*r0.*d2rdulduk;
            d2Ddulduk = 2*R^2*(dqdul.*r0.*r.*f - ...
              q.*(2*r0.*f.*drdul + 2.*r.*f.*dr0dul + r0.*r.*dfdul/2))./denom0; 
            H(l,k) = H(l,k) + n/(2*R)*s*(16*R^4*(dDdul.*dDduk)./denom1 - ...
              24*R^4*(dDdul.*dDduk)./denom2 - 8*R^4*(d2Ddulduk)./denom3);
            H(k,l) = H(l,k);
          end
        end
      end
    end
  end
  
% AK contribution:
  g1 = bsxfun(@times,kz'+w,cos(t*kk'));
  g2 = bsxfun(@times,kz'+w,sin(t*kk'));
  g = (g1*u - g2*v).^2 + (g2*u + g1*v).^2;
  h = (R^2 + (a0*u + b0*v).^2 + (-b0*u + a0*v).^2).^2;
  h3 = h.^3;
  dgdu = 2*(bsxfun(@times,g1*u - g2*v,g1) + bsxfun(@times,g2*u + g1*v,g2));
  dhdu = bsxfun(@times,a0*u + b0*v,a0) + bsxfun(@times,-b0*u + a0*v,-b0);
  dhdu = bsxfun(@times,4*sqrt(h),dhdu);
  G(1:N) = G(1:N) + 4*R^4*(n/2)*(bsxfun(@rdivide,dgdu,h)'*s' - ...
    bsxfun(@rdivide,dhdu,h.^2./g)'*s');
  dgdv = 2*(bsxfun(@times,g2*u + g1*v,g1) + bsxfun(@times,g1*u - g2*v,-g2));
  dhdv = bsxfun(@times,a0*u + b0*v,b0) + bsxfun(@times,-b0*u + a0*v,a0);
  dhdv = bsxfun(@times,4*sqrt(h),dhdv);
  G(N+1:2*N) = G(N+1:2*N) + 4*R^4*(n/2)*(bsxfun(@rdivide,dgdv,h)'*s' -...
     bsxfun(@rdivide,dhdv,h.^2./g)'*s');
  if nargout > 1
    for l = 1:N 
      dgdul = dgdu(:,l);
      dgdvl = dgdv(:,l);
      dhdul = dhdu(:,l);
      dhdvl = dhdv(:,l);
      for k = 1:N 
        dgduk = dgdu(:,k);
        dgdvk = dgdv(:,k);
        dhduk = dhdu(:,k);
        dhdvk = dhdv(:,k);
        p = h.*dgdvk - g.*dhdvk;
        d2gduldvk = 2*(kz(k)+w)*(kz(l)+w)*sin((kk(l)-kk(k))*t);
        d2hduldvk = 1./(2*h).*dhdul.*dhdvk + 4*sqrt(h).*sin((kk(l)-kk(k))*t);
        dpdul = dhdul.*dgdvk + h.*d2gduldvk - dgdul.*dhdvk - g.*d2hduldvk;
        H(l,k+N) = H(l,k+N) + (2*R^4*n)*s*((dpdul.*h - 2*p.*dhdul)./h3);
        if k <= l
          d2gdvldvk = 2*(kz(k)+w)*(kz(l)+w)*cos((kk(k)-kk(l))*t);
          d2hdvldvk = 1./(2*h).*dhdvl.*dhdvk+4*sqrt(h).*cos((kk(k)-kk(l))*t);
          dpdvl = dhdvl.*dgdvk + h.*d2gdvldvk - dgdvl.*dhdvk - g.*d2hdvldvk;
          H(l+N,k+N) = H(l+N,k+N) + (2*R^4*n)*s*((dpdvl.*h - 2*p.*dhdvl)./h3);
          H(k+N,l+N) = H(l+N,k+N);
          p = h.*dgduk - g.*dhduk;
          d2hdulduk = 1./(2*h).*dhdul.*dhduk+4*sqrt(h).*cos((kk(k)-kk(l))*t);
          dpdul = dhdul.*dgduk + h.*d2gdvldvk - dgdul.*dhduk - g.*d2hdulduk;
          H(l,k) = H(l,k) + (2*R^4*n)*s*((dpdul.*h - 2*p.*dhdul)./h3);
          H(k,l) = H(l,k);
        end
      end
    end
  end
  
% Get rid of the constant term and use symmetry:
  if ( mod(N, 2) == 0 )
    idx1 = N/2;
    idx2 = N+N/2;
  else
    idx1 = (N-1)/2;
    idx2 = N+(N-1)/2;
  end
  G = [G(1:idx1); G(idx1+2:idx2); G(idx2+2:end)];
  if nargout > 1
    H = [H(1:idx1,1:idx1), H(1:idx1,idx1+2:idx2), H(1:idx1,idx2+2:end);...
       H(idx1+2:idx2,1:idx1), H(idx1+2:idx2,idx1+2:idx2), ...
       H(idx1+2:idx2,idx2+2:end); H(idx2+2:end,1:idx1), ...
       H(idx2+2:end,idx1+2:idx2), H(idx2+2:end,idx2+2:end)];
    H(N:end,1:N-1) = H(1:(N-1),N:end)';
  end
  
end