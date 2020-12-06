function plotonplane(q,n,w)

vscale = max(max(real(q)),max(imag(q)));
xmax = 1.2*vscale; ymax = xmax;
LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
lw = 2; ms = 30; fs = 18;
dom = [0 2*pi]; T = 2*pi; dt = .1;
xStars = 2*xmax*rand(250,1)-xmax;
yStars = 2*ymax*rand(250,1)-ymax;
for t = dt:dt:10*T
    clf, hold off
    fill(xmax*[-1 1 1 -1 -1],ymax*[-1 -1 1 1 -1],'k')
    axis equal, axis([-xmax xmax -ymax ymax]), box on
    set(gca,FS,fs)
    hold on, plot(xStars,yStars,'.w',MS,4)
    if w ~= 0
        q = chebfun(@(t)exp(1i*w*dt).*q(t),dom,length(q),'trig');
    end
    if t < T
        hold on, plot(q,'r',LW,lw)
    end
    for j = 1:n
        hold on, plot(q(t+2*pi*j/n),'.',MS,ms)
    end
    drawnow
end

end