function plotonsphere(q,n,w,R)

if nargin < 3
    w = 0;
    R = 1;
end
LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
lw = 2; ms = 30; fs = 18;
XTick = [-R,0,R]; YTick = XTick; ZTick = XTick;
N = 500;
Ns = 100;
ns = 50;
[Xs,Ys,Zs] = sphere(Ns);
Xs = R*Xs; Ys = R*Ys; Zs = R*Zs;
XXs = Xs(1:ns:end,:); YYs = Ys(1:ns:end,:); ZZs = Zs(1:ns:end,:);
XXs = XXs'; YYs = YYs'; ZZs = ZZs';
T = 2*pi;
dt = .1;
tt = trigpts(N, [0 2*pi]);
z = q(tt);
[X,Y,Z] = plane2sphere(z,R);
for t = dt:dt:10*T
  clf, hold off
  surf(Xs,Ys,Zs,'EdgeColor','none'), alpha(.15)
  hold on, plot3(Xs(:,1:ns:end),Ys(:,1:ns:end),Zs(:,1:ns:end),'k',LW,1e-10)
  hold on, plot3(XXs,YYs,ZZs,'k',LW,1e-10), view(-45.5,46)
  if t > 5*T
    view(0,90)
  end
  if t < 3*T
    if w ~= 0
      M = [X';Y';Z'];
      Rot = [cos(w*dt),-sin(w*dt),0;sin(w*dt),cos(w*dt),0;0,0,1];
      M = Rot*M;
      X = M(1,:)'; Y = M(2,:)'; Z = M(3,:)';
    end
  hold on, plot3(X,Y,Z,'r',LW,lw)
  end
  for j=1:n
    z = q(t+2*pi*j/n);
    [XX,YY,ZZ] = plane2sphere(z,R);
    if w ~= 0
      M = [XX';YY';ZZ'];
      Rot = [cos(w*t),-sin(w*t),0;sin(w*t),cos(w*t),0;0,0,1];
      M = Rot*M;
      XX = M(1,:)'; YY = M(2,:)'; ZZ = M(3,:)';
    end
    hold on, plot3(XX,YY,ZZ,'.',MS,ms)
  end
  box on, set(gca,FS,fs), set(gca,'XTick',XTick,'YTick',YTick,'ZTick',ZTick)
  axis([-R R -R R -R R]), axis equal
  drawnow
end

end

function [x1,x2,x3] = plane2sphere(z,R)
  if nargin < 2
    R = 1;
  end
  x1 = 2*R^2*real(z)./(R^2 + abs(z).^2);  
  x2 = 2*R^2*imag(z)./(R^2 + abs(z).^2);
  x3 = R*(-R^2 + abs(z).^2)./(R^2 + abs(z).^2); 
end
