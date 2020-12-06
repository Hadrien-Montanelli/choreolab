% choreo.m - find choreographies on the plane via minimization of the action and 
%           using a hand-drawn initial guess.

close all
LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
lw = 2; ms = 30; fs = 18;
format long, format compact
n = input('How many bodies? (e.g. 5) '); 
w = input('Angular velocity? (non-integer, e.g. 1.2) ');
k = 15;
N = k*n;
dom = [0 2*pi];

while 1   

% Hand-drawn initial guess:
  xmax = 2.5; ymax = 2.5;
  clf, hold off, axis equal, axis([-xmax xmax -ymax ymax]), box on
  set(gca,FS,fs), title('Draw a curve!')
  h = imfreehand;
  z = getPosition(h);
  delete(h)
  z = z(:,1) + 1i*z(:,2);
  q0 = chebfun(z,dom,'trig');
  q0 = chebfun(q0,dom,N,'trig');
  c0 = trigcoeffs(q0);
  c0(1+floor(N/2)) = 0;
  q0 = chebfun(c0,dom,'coeffs','trig');
    
% Solve the problem:
  hold on, plot(q0,'.-b',LW,lw), drawnow
  c0 = trigcoeffs(q0);
  c0 = [real(c0);imag(c0)];
  [A0,G0] = actiongradeval(c0,n,w);
  fprintf('\nInitial acion: %.6f\n',A0)
  options = optimoptions('fminunc');
  options.Algorithm = 'quasi-newton';
  options.HessUpdate = 'bfgs'; 
  options.GradObj = 'on';
  options.Display = 'off';
  [c,A,~, ~,G] = fminunc(@(x)actiongradeval(x,n,w),c0,options);
  fprintf('Action after optimization: %.6f\n',A)
  fprintf('Norm of the gradient: %.3e\n',norm(G)/norm(G0))

% Plot the result:
  c = c(1:N) + 1i*c(N+1:2*N);
  c(1+floor(N/2)) = 0;
  q = chebfun(c,dom,'coeffs','trig');
  hold on, plot(q,'r',LW,lw)
  hold on, plot(q(2*pi*(0:n-1)/n),'.k',MS,ms), title('')
  s = input('Want to see the planets dancing? (Yes=1, No=0) '); 
  if s == 1
    hold off, plotonplane(q,n,w), pause
  end
  
end
