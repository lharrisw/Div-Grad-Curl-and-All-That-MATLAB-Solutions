%% I-1 and I-6
%   Plot the following vector fields:
%   
%   (a) iy + jx           (e) jx
%   (b) (i+j)/sqrt(2)     (f) (iy + jx)/sqrt(x^2 + y^2),(x,y) =/= (0,0)
%   (c) ix - jy           (g) iy +jxy
%   (d) iy                (h) i + jy

clear
clc
close all

[x,y] = meshgrid(-3:0.5:3);
xx = linspace(-3,3);

% part a
u = y;
v = x;

% I-6
yy_p = @(c)  sqrt(xx.^2-c);
yy_m = @(c) -sqrt(xx.^2-c);

figure
hold on
quiver(x,y,u,v)
for j = -3:3
    yp = yy_p(j);
    ym = yy_m(j);
    plot(real(xx(yp>0)),real(yp(yp>0)),'b.')
    plot(real(xx(ym<0)),real(ym(ym<0)),'b.')
end


% % part b
% u = ones(size(x))/sqrt(2);
% v = u;
% 
% % I-6
% yy = @(c) xx + c;
% 
% figure
% hold on
% quiver(x,y,u,v)
% for j = -3:3
%     yc = yy(j);
%     plot(xx,yc,'b.')
% end
% xlim([-3 3])
% ylim([-3 3])

% % part c
% u = x;
% v = -y;
% 
% % I-6
% yy = @(c) c./xx;
% 
% figure
% hold on
% quiver(x,y,u,v)
% for j = -3:3
%     if j == 0
%         yc = zeros(size(xx));
%         plot(xx,yc,'b.')
%     end
%     yc = yy(j);
%     plot(xx,yc,'b.')
% end
% plot(0*xx,xx,'b.')
% xlim([-3 3])
% ylim([-3 3])

% % part d
% u = y;
% v = 0*u;
% 
% % I-6
% yy = @(c) c;
% 
% figure
% hold on
% quiver(x,y,u,v)
% for j = -3:3
%     yc = yy(j);
%     plot(xx,yc,'b.')
% end

% % part e
% u = zeros(size(x));
% v = x;
% 
% % I-6
% figure
% hold on
% quiver(x,y,u,v)
% for j = -3:3
%     plot(j*ones(size(xx)),xx,'b.')
% end

% % part f
% u = y./sqrt(x.^2+y.^2);
% v = x./sqrt(x.^2+y.^2);
% 
% % I-6
% yy_p = @(c)  sqrt(xx.^2-c);
% yy_m = @(c) -sqrt(xx.^2-c);
% 
% figure
% hold on
% quiver(x,y,u,v)
% for j = -3:3
%     yp = yy_p(j);
%     ym = yy_m(j);
%     plot(real(xx(yp>0)),real(yp(yp>0)),'b.')
%     plot(real(xx(ym<0)),real(ym(ym<0)),'b.')
% end

% % part g
% u = y;
% v = x.*y;
% 
% % I-6
% yy = @(c) xx.^2/2 + c;
% 
% figure
% hold on
% quiver(x,y,u,v)
% for j = -3:3
%     yc = yy(j);
%     plot(xx,yc,'b.')
% end
% xlim([-3 3])
% ylim([-3 3])

% % part h
% u = ones(size(x));
% v = y;
% 
% % I-6
% yy = @(c) c*exp(xx);
% 
% figure
% hold on
% quiver(x,y,u,v)
% for j = -3:3
%     yc = yy(j);
%     plot(xx,yc,'b.')
% end
% xlim([-3 3])
% ylim([-3 3])

%% I-2
%   Sketch the electric field of a unit positive charge situated at the
%   origin.

clear
clc
close all

% R2
[x,y] = meshgrid(linspace(-3,3,25));
u = x./(x.^2+y.^2).^1.5;
v = y./(x.^2+y.^2).^1.5;

figure
hold on
quiver(x,y,u,v)

% R3
[x,y,z] = meshgrid(linspace(-3,3,25));
u = x./(x.^2+y.^2+z.^2).^1.5;
v = y./(x.^2+y.^2+z.^2).^1.5;
w = z./(x.^2+y.^2+z.^2).^1.5;

figure
hold on
quiver3(x,y,z,u,v,w)
view(3)

%% I-3
%   (a) Write a formula for a vector function in two dimensions, which is
%   in the positive radial direction and whose magnitude is 1.
%   (b) Write a formula for a vector function in two dimensions whose
%   direction makes an angle of 45 degrees with the x-axis and whose
%   magnitude at any point (x,y) is (x+y)^2.
%   (c) Write a formula for a vector function in two dimensions whose
%   direction is tangential and whose magnitude at any point (x,y) is equal
%   to its distance from the origin.
%   (d) Write a formula for a vector function in three dimensions which is
%   in the positive radial direction and whose magnitude is 1.

clear
clc
close all

% grid
% [x,y] = meshgrid(linspace(-3,3,25));

% % part a
% u = x./sqrt(x.^2+y.^2);
% v = y./sqrt(x.^2+y.^2);

% % part b
% u = sqrt(2)*(x+y).^2/2;
% v = sqrt(2)*(x+y).^2/2;

% % part c
% u = -y;
% v = x;

% figure
% hold on
% quiver(x,y,u,v)

% part d
[x,y,z] = meshgrid(linspace(-3,3,10));
u = x./sqrt(x.^2+y.^2+z.^2);
v = y./sqrt(x.^2+y.^2+z.^2);
w = z./sqrt(x.^2+y.^2+z.^2);

figure
hold on
quiver3(x,y,z,u,v,w)
view(3)

%% I-4
%   An object moves in the xy-plane in such a way that its position vector
%   r is given as a function of time t by r = i(a cos(wt)) + j(b sin(wt)),
%   where a, b, and w are constants.
%   (a) how far is the object from the origin at any time t?
%   (b) Find the object's velocity and acceleration as functions of time.
%   (c) Show that the obhect moves on the elliptical path
%   (x/a)^2+(y/b)^2=1.

clear
clc
close all

t = linspace(0,1);
a = 3; % x-coordinate of the center
b = 5; % y-coordinate of the center
w = 2*pi; % frequency

% position
x = a*cos(w*t); % parameterization of x
y = b*sin(w*t); % parameterization of y

% velocity
xdot = -a*w*sin(w*t); % x-component of velocity
ydot = b*w*cos(w*t); % y-component of velocity

% acceleration
xddot = -a*w^2*cos(w*t); % x-component of acceleration
yddot = -b*w^2*sin(w*t); % y-component of acceleration

figure
hold on
plot(x,y,'b') % plot of the parametric curve
q1 = quiver(x,y,xdot,ydot,0.075,'b-','linewidth',1.5); % plot of the velocity vectors
q2 = quiver(x,y,xddot,yddot,0.01,'r-','linewidth',1.5); % plot of the acceleration vectors
% quiver(zeros(size(x)),zeros(size(y)),x,y,0) % plot of the position vectors

for k = 1:length(t)
    % reassignments for position
    xx = x(k);
    yy = y(k);
    plot(xx,yy,'b.','markersize',10);
    % reassignments for velocity vectors
    xxdot = xdot(k);
    yydot = ydot(k);
    % reassignments for acceleration vectors
    xxddot = xddot(k);
    yyddot = yddot(k);
    % updates for velocity
    q1.XData = xx;
    q1.YData = yy;
    q1.UData = xxdot;
    q1.VData = yydot;
    % updates for acceleration
    q2.XData = xx;
    q2.YData = yy;
    q2.UData = xxddot;
    q2.VData = yyddot;
    xlim([min(x)-1 max(x)+1])
    ylim([min(y)-1 max(y)+1])
    hold on
    grid on
    pause(0.1)
end

%% I-5
%   A charge +1 is situated at the point (1,0,0) and a charge -1 is
%   situated at the point (-1,0,0). Find the electric field of these two
%   charges at an arbitrary point (0,y,0) on the y-axis.

clear
clc
close all

% grid
[x,y,z] = meshgrid(linspace(-3,3,10));

% electric field at (0,k,0) due to +1 charge at (1,0,0)
E_px = @(k) -ones(size(x))./(1+ones(size(x))*k^2).^1.5; % x-component
E_py = @(k) ones(size(x))*k./(1+ones(size(x))*k^2).^1.5; % y-component
E_pz = @(k)  zeros(size(x)); % z-component

% electric field at (0,k,0) dut to -1 charge at (-1,0,0)
E_mx = @(k)  E_px(k); % x-component
E_my = @(k) -E_py(k); % y-component
E_mz = @(k)  E_pz(k); % z-component

% total field at (0,k,0)
Ex = @(k) E_px(k) + E_mx(k);
Ey = @(k) E_py(k) + E_my(k);
Ez = @(k) E_pz(k) + E_mz(k);

k0 = 1;

figure
hold on
grid on
plot3(1,0,0,'k.','markersize',30)
plot3(-1,0,0,'r.','markersize',30)
plot3(0,k0,0,'b.','markersize',30)
quiver3(0*x,k0*ones(size(y)),0*z,E_px(k0),E_py(k0),E_pz(k0),0,'k')
quiver3(0*x,k0*ones(size(y)),0*z,E_mx(k0),E_my(k0),E_mz(k0),0,'r')
quiver3(0*x,k0*ones(size(y)),0*z,Ex(k0),Ey(k0),Ez(k0),0,'b')
xlabel('x')
ylabel('y')
zlabel('z')
view(3)