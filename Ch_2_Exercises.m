%% II-1
%   Find the unit vector n normal to each of the following surfaces
%   (a) z = 2 - x - y
%   (b) z = (x^2+y^2)^1/2
%   (c) z = (1-x^2)^1/2
%   (d) z = x^2+y^2
%   (e) z = (1-x^2/a^2-y^2/a^2)^1/2

clear
clc
close all

[x,y] = meshgrid(linspace(-3,3,40));

xx = 1.6;
yy = 1;

% % (a)
% % function
% z = 2-x-y;

% % position vector
% R = @(x,y) [x;y;2-x-y];
% r = R(xx,yy);

% %  normal vector for a single point
% n = [1 1 1]/sqrt(3);

% % normal vectors over surface
% u = 1/sqrt(3)*ones(size(x));
% v = u;
% w = v;

% % (b)
% % function
% z = sqrt(x.^2+y.^2);

% % position vector
% R = @(x,y) [x;y;sqrt(x.^2+y.^2)];
% r = R(xx,yy);

% % normal vector for a point
% N = @(x,y) -[-x;-y;sqrt(x.^2+y.^2)]./(sqrt(2)*sqrt(x.^2+y.^2));
% n = N(xx,yy);

% % normal vectors for the surface
% u = x./(sqrt(2)*sqrt(x.^2+y.^2));
% v = y./(sqrt(2)*sqrt(x.^2+y.^2));
% w = -z./(sqrt(2)*sqrt(x.^2+y.^2));

% % (c)
% % function
% z = real(sqrt(1-x.^2));

% % position vector
% R = @(x,y) [x;y;sqrt(1-x.^2)];
% r = R(xx,yy);

% % normal vector for single point
% N = @(x,y) [x;0;sqrt(1-x.^2)];
% n = N(xx,yy);

% % normal vectors over surface
% u = x;
% v = 0*u;
% w = z;

% % (d)
% %  function
% z = x.^2+y.^2;

% % position vector
% R = @(x,y) [x;y;x.^2+y.^2];
% r = R(xx,yy);

% % normal vector for single point
% N = @(x,y) [-2*x;-2*y;1]./sqrt(4*(x.^2+y.^2)+1);
% n = N(xx,yy);

% % normal vectors over surface
% u = 2*x./sqrt(4*(x.^2+y.^2)+1);
% v = 2*y./sqrt(4*(x.^2+y.^2)+1);
% w = -1./sqrt(4*(x.^2+y.^2)+1);

% % (e)
% a = 2; % radius

% % function
% z = real(sqrt(1-x.^2/a^2-y.^2/a^2));

% % position vector
% R = @(x,y) [x;y;real(sqrt(1-x.^2/a^2-y.^2/a^2))];
% r = R(xx,yy);

% % normal vector for a single point
% N = @(x,y) [x;y;a^2*real(sqrt(1-x.^2/a^2-y.^2/a^2))]/a...
%     ./sqrt(1+(a^2-1)*real(sqrt(1-x.^2/a^2-y.^2/a^2)));
% n = N(xx,yy);

% % normal vectors over the surface
% u = x/a./sqrt(1+(a^2-1)*real(sqrt(1-x.^2/a^2-y.^2/a^2)));
% v = y/a./sqrt(1+(a^2-1)*real(sqrt(1-x.^2/a^2-y.^2/a^2)));
% w = a*z./sqrt(1+(a^2-1)*real(sqrt(1-x.^2/a^2-y.^2/a^2)));

figure
hold on
grid on
mesh(x,y,z)
% quiver3(r(1),r(2),r(3),n(1),n(2),n(3),'k-','linewidth',2)
quiver3(x,y,z,u,v,w)
axis equal
view(3)

%% II-2
%   Show that the unit vector normal to the plane ax + by + cz = d is given
%   by n = +/-(ai + bj + ck)/sqrt(a^2+b^2+c^2).

clear
clc
close all

[x,y] = meshgrid(linspace(-3,3,25));
% components of normal vector
a = 1;
b = -1;
c = 1;

% length of normal vector
L = norm([a b c],2);
d = 0;
if c == 0
    [x,z] = meshgrid(linspace(-3,3,25));
    % equation of a plane
    y = d/b - a*x/b; 
    % components of unit normal
    u = a*ones(size(x))/L;
    v = b*ones(size(x))/L;
    w = c*ones(size(x))/L;
    
    figure
    hold on
    grid on
    mesh(x,y,z)
    quiver3(x,y,z,u,v,w)
    axis equal
    view(3)
else
    % equation of a plane
    z = d/c-a*x/c-b*y/c;    
    % components of unit normal
    u = a*ones(size(x))/L;
    v = b*ones(size(x))/L;
    w = c*ones(size(x))/L;

    figure
    hold on
    grid on
    mesh(x,y,z)
    quiver3(x,y,z,u,v,w)
    axis equal
    view(3)
end

%% II-4
clear
clc
close all

% % (a)
% [x,z] = meshgrid(linspace(-3,3,25));

% % function of integration and surface
% y = 1-x-z;
% 
% % unit normal
% nx = ones(size(x))/sqrt(3);
% ny = ones(size(y))/sqrt(3);
% nz = ones(size(x))/sqrt(3);
% 
% figure
% hold on
% grid on
% mesh(x,y,z)
% quiver3(x,y,z,nx,ny,nz,0)
% axis equal
% view(3)

% % (b)
% % surface 1: paraboloid
% [xx,yy] = meshgrid(linspace(-1,1,25));
% zz = xx.^2+yy.^2;
% 
% % function of integration
% G = 1./sqrt(4*(xx.^2+yy.^2)+1);
% 
% % surface 2: circle
% [r,theta] = meshgrid(linspace(0,1,25),linspace(0,2*pi,25));
% x = r.*cos(theta);
% y = r.*sin(theta);
% z = ones(size(x));
% 
% % normal vector to the surface
% nx = -2*xx./sqrt(4*(xx.^2+yy.^2)+1);
% ny = -2*yy./sqrt(4*(xx.^2+yy.^2)+1);
% nz = 1./sqrt(4*(x.^2+y.^2)+1);
% 
% figure
% hold on
% grid on
% mesh(x,y,z)
% mesh(xx,yy,zz)
% % mesh(xx,yy,G)
% quiver3(xx,yy,zz,nx,ny,nz)
% zlim([0 1])
% axis equal
% view(3)

% % (c)
% [x,y] = meshgrid(linspace(-3,3,30));
% 
% % function of integration
% G = real((1-x.^2-y.^2).^1.5);
% 
% % surface
% z = real((1-x.^2-y.^2).^0.5);
% 
% % normal vector
% nx = x;
% ny = y;
% nz = z;
% 
% figure
% hold on
% grid on
% mesh(x,y,z)
% mesh(x,y,G)
% quiver3(x,y,z,nx,ny,nz)
% axis equal
% view(3)

%% II-5

clear
clc
close all

% % (a)
% [x,z] = meshgrid(linspace(-3,3,30));
% 
% % vector field
% Fx = x;
% Fy = 0*x;
% Fz = -z;
% 
% % surface
% y = 2-x-2*z;
% 
% % normal vectors to the surface
% nx = ones(size(x))/sqrt(6);
% ny = nx;
% nz = 2*ny;
% 
% % normal component of field
% Fnx = nx.*(Fx.*nx+Fy.*ny+Fz.*nz);
% Fny = ny.*(Fx.*nx+Fy.*ny+Fz.*nz);
% Fnz = nz.*(Fx.*nx+Fy.*ny+Fz.*nz);
% 
% figure
% hold on
% grid on
% mesh(x,y,z)
% quiver3(x,y,z,Fx,Fy,Fz,0)
% quiver3(x,y,z,nx,ny,nz,0)
% quiver3(x,y,z,Fnx,Fny,Fnz,0)
% axis equal
% view(3)

% % (b)
% [x,y] = meshgrid(linspace(-3,3,25));
% 
% % radius and surface
% a = 2;
% z = real(sqrt(a^2-x.^2-y.^2));
% 
% % vector field
% Fx = x;
% Fy = y;
% Fz = z;
% 
% % normal vectors to the surface
% nx = x/a;
% ny = y/a;
% nz = z/a;
% 
% % normal component of field
% Fnx = nx.*(Fx.*nx+Fy.*ny+Fz.*nz);
% Fny = ny.*(Fx.*nx+Fy.*ny+Fz.*nz);
% Fnz = nz.*(Fx.*nx+Fy.*ny+Fz.*nz);
% 
% figure
% hold on
% grid on
% mesh(x,y,z)
% quiver3(x,y,z,Fx,Fy,Fz,0)
% quiver3(x,y,z,nx,ny,nz,0)
% quiver3(x,y,z,Fnx,Fny,Fnz,0)
% axis equal
% xlim([-3 3])
% ylim([-3 3])
% view(3)

% % (c)
% [x,y] = meshgrid(linspace(-1,1,25));
% z = 1-x.^2-y.^2;
% 
% nx = 2*x./sqrt(4*(x.^2+y.^2)+1);
% ny = 2*y./sqrt(4*(x.^2+y.^2)+1);
% nz = 1./sqrt(4*(x.^2+y.^2)+1);
% 
% figure
% hold on
% grid on
% mesh(x,y,z)
% quiver3(x,y,z,nx,ny,nz)
% zlim([0 1.5])
% axis equal
% view(3)

%% II-6

clear
clc
close all

R = 1;
s0 = 1;
[x,y] = meshgrid(linspace(-R,R,50));

% surface
z = real(sqrt(R^2-x.^2-y.^2));

% mass distribution
s = s0*(x.^2+y.^2)/R^2;

figure
hold on
grid on
mesh(x,y,z)
% mesh(x,y,s)
axis equal
view(3)

%% II-9

clear
clc
close all

n = 30; % grid size
h = 3; % height
r = 2; % radius
l = 1; % lambda
[x,y] = meshgrid(linspace(-h/2,h/2,n),linspace(-r,r,n));
[yy,zz] = meshgrid(linspace(-r,r,n),linspace(0,r,n));

% surface
z = real(sqrt(r^2-y.^2));
xx = ones(size(x))*h/2;

% electric field
Ex = l*x;
Ey = l*y;
Ez = l*zeros(size(y));

figure
hold on
grid on
mesh(x,y,z)
mesh(xx,yy,zz)
mesh(-xx,yy,zz)
quiver3(x,y,z,Ex,Ey,Ez)
view(3)

%% II-10

clear
clc
close all

R = 2;
h = 5;
[x,z] = meshgrid(linspace(-R,R,25),linspace(0,h,25));

% surface
y = real(sqrt(R^2-x.^2));

% normal vectors to surface
nx = x./R;
ny = y/R;
nz = 0*ny;

figure
hold on
grid on
mesh(x,y,z)
mesh(x,-y,z)
quiver3(x,y,z,nx,ny,nz)
quiver3(x,-y,z,nx,-ny,nz)
axis equal
view(3)

%% II-14
clear
clc
close all

[x,y,z] = meshgrid(linspace(-2,2,10));
% a

u = x.^2;
v = y.^2;
w = z.^2;

% divergence
div = 2*(x+y+z);

figure
hold on
grid on
quiver3(x,y,z,u,v,w)
view(3)

figure
slice(x,y,z,div,[-2,0,2],0,0)

%% II-24
clear
clc
close all

[x,y] = meshgrid(linspace(-2,2,25));
h = 1;
z = h-sqrt(x.^2+y.^2);
u = x/sqrt(2)./sqrt(x.^2+y.^2);
v = y/sqrt(2)./sqrt(x.^2+y.^2);
w = ones(size(z))/sqrt(2);

figure
hold on
grid on
mesh(x,y,z)
% quiver3(x*0,y*0,z*0,u,v,w)
quiver3(x,y,z,u,v,w)
zlim([0 h])
axis equal
view(3)

%% II-26
clear
clc
close all

h = 1;
R = 1;
L = sqrt(R^2+h^2);
[r,t] = meshgrid(linspace(0,R,25),linspace(0,2*pi,25));
x = r.*cos(t);
y = r.*sin(t);
z = h*sqrt(x.^2+y.^2)/R;
u = x*h./sqrt(x.^2+y.^2)/L;
v = y*h./sqrt(x.^2+y.^2)/L;
w = -ones(size(z))*R/L;

figure
hold on
grid on
mesh(x,y,z)
quiver3(x,y,z,u,v,w)
axis equal
view(3)