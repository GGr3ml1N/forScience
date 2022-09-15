clear;
close all;
clc;

%% Инициализация входных параметров
p = 50;
a = -1;
b = 1;
R1 = linspace(a,b,p);
[X,Y] = meshgrid(R1,R1);
[theta1,rho] = cart2pol(X,Y);
NA = 0.95;
lmd = 0.532;
k = 2*pi/lmd;
f = 1;
t = 0;
B = 1;
R2 = 1;
alpha = asin(NA);
% graph_z = 0:0.05:2;


Ex = zeros(p);
Ey = zeros(p);
Ez = zeros(p);
Hx = zeros(p);
Hy = zeros(p);
Hz = zeros(p);
Ex_1D = zeros(size(R1));
Ey_1D = zeros(size(R1));
Ez_1D = zeros(size(R1));
Hx_1D = zeros(size(R1));
Hy_1D = zeros(size(R1));
Hz_1D = zeros(size(R1));
X1 = zeros(1,41);
Y1 = zeros(1,41);


% % % %  X1 = [-2.212182563139366e-10, -0.0028, -0.0055, -0.0079, -0.0098, -0.0111, ];
% % % %  Y1 = [4.456445262045631e-09,   0.0410,  0.0815,  0.1213,  0.1599,  0.1970,  ];
i = 1;
% Расчеты
tic
for z = (0:0.05:2)*lmd

parfor n = 1:p
    for m = 1:p
    Ex(n,m) = -1i.*f/(lmd).*integral2(@(theta,phi) (sqrt(cos(theta))).*((sin(phi).^2).*(cos(2.*phi))+(cos(phi).^2).*(cos(theta)).*(cos(2.*phi))+(sin(phi)).*(cos(phi)).*(cos(theta)).*(sin(2.*phi))-(sin(phi)).*(cos(phi)).*(sin(2.*phi))).*exp(1i.*k.*(rho(n,m).*sin(theta).*cos(phi-theta1(n,m))+z.*cos(theta))).*sin(theta),0,alpha,0,2*pi);
    end
end

parfor n = 1:p
    for m = 1:p
    Ey(n,m) = -1i.*f/(lmd).*integral2(@(theta,phi) (sqrt(cos(theta))).*(sin(phi).*cos(theta).*cos(phi)-sin(phi).*cos(phi).*cos(2.*phi)+(cos(phi).^2).*sin(2.*phi)).*exp(1i.*k.*(rho(n,m).*sin(theta).*cos(phi-theta1(n,m))+z.*cos(theta))).*sin(theta),0,alpha,0,2*pi);
    end
end

parfor n = 1:p
    for m = 1:p
    Ez(n,m) = -1i.*f/(lmd).*integral2(@(theta,phi) (sqrt(cos(theta))).*(-sin(theta).*cos(phi)).*exp(1i.*k.*(rho(n,m).*sin(theta).*cos(phi-theta1(n,m))+z.*cos(theta))).*sin(theta),0,alpha,0,2*pi);
    end
end

parfor n = 1:p
    for m = 1:p
    Hx(n,m) = -1i.*f/(lmd).*integral2(@(theta,phi) (sqrt(cos(theta))).*(-sin(phi).^2.*sin(2.*phi)-(cos(phi).^2).*cos(theta).*sin(2.*phi)+sin(phi).*cos(phi).*cos(2.*phi).*cos(theta)-sin(phi).*cos(phi).*cos(2.*phi)).*exp(1i.*k.*(rho(n,m).*sin(theta).*cos(phi-theta1(n,m))+z.*cos(theta))).*sin(theta),0,alpha,0,2*pi);
    end
end

parfor n = 1:p
    for m = 1:p
    Hy(n,m) = -1i.*f/(lmd).*integral2(@(theta,phi) (sqrt(cos(theta))).*(-sin(phi).*cos(phi).*sin(2.*phi).*cos(theta)+sin(phi).*cos(phi).*sin(2.*phi)+(cos(phi).^2).*cos(2.*phi)+(sin(phi).^2).*cos(theta).*cos(2.*phi)).*exp(1i.*k.*(rho(n,m).*sin(theta).*cos(phi-theta1(n,m))+z.*cos(theta))).*sin(theta),0,alpha,0,2*pi);
    end
end

parfor n = 1:p
    for m = 1:p
    Hz(n,m) = -1i.*f/(lmd).*integral2(@(theta,phi) (sqrt(cos(theta))).*(sin(theta).*sin(phi)).*exp(1i.*k.*(rho(n,m).*sin(theta).*cos(phi-theta1(n,m))+z.*cos(theta))).*sin(theta),0,alpha,0,2*pi);
    end
end

S = abs((Ey.*conj(Hz)-Ez.*conj(Hy))-(Ex.*conj(Hz)-Ez.*conj(Hx))+(Ex.*conj(Hy)-Ey.*conj(Hx)));
Sx = (1/2) * real(Ey.*conj(Hz)-Ez.*conj(Hy));
Sy = (1/2) * real(Ex.*conj(Hz)-Ez.*conj(Hx));
Sz = (1/2) * real(Ex.*conj(Hy)-Ey.*conj(Hx));

x = min(Sx(p/2,p/2:end));
y = max(Sx(p/2,p/2:end));


 X1(i)=x;
 Y1(i)=y;
 i=i+1;
end
time = toc;
% Графическое представление
figure('Name','Ex','NumberTitle','off'); 
imagesc(R1, R1, Ex.*conj(Ex)); 
xlabel('x, \mum'); 
ylabel('y, \mum'); 
axis square xy; 
colorbar

figure('Name','Ey','NumberTitle','off'); 
imagesc(R1, R1, Ey.*conj(Ey)); 
xlabel('x, \mum'); 
ylabel('y, \mum'); 
axis square xy; 
colorbar

figure('Name','Ez','NumberTitle','off'); 
imagesc(R1, R1, Ez.*conj(Ez)); 
xlabel('x, \mum'); 
ylabel('y, \mum'); 
axis square xy; 
colorbar

figure('Name','E1','NumberTitle','off');
imagesc(R1, R1, (Ex.*conj(Ex) + Ey.*conj(Ey) + Ez.*conj(Ez)));
xlabel('x, \mum');
ylabel('y, \mum');
axis square xy;
colorbar

figure('Name','Hx','NumberTitle','off'); 
imagesc(R1, R1, Hx.*conj(Hx)); 
xlabel('x, \mum'); 
ylabel('y, \mum'); 
axis square xy; 
colorbar

figure('Name','Hy','NumberTitle','off'); 
imagesc(R1, R1, Hy.*conj(Hy)); 
xlabel('x, \mum'); 
ylabel('y, \mum'); 
axis square xy; 
colorbar

figure('Name','Hz','NumberTitle','off'); 
imagesc(R1, R1, Hz.*conj(Hz)); 
xlabel('x, \mum'); 
ylabel('y, \mum'); 
axis square xy; 
colorbar

figure('Name','H1','NumberTitle','off');
imagesc(R1, R1, (Hx.*conj(Hx) + Hy.*conj(Hy) + Hz.*conj(Hz)));
xlabel('x, \mum');
ylabel('y, \mum');
axis square xy;
colorbar

figure('Name','Ex_1D','NumberTitle','off'); 
plot(R1, Ex(p/2,:).*conj(Ex(p/2,:)));
xlabel('x, \mum'); 
ylabel('Ez, V/m'); 
axis square xy; 

figure('Name','Ey_1D','NumberTitle','off'); 
plot(R1, Ey(p/2,:).*conj(Ey(p/2,:)));
xlabel('x, \mum'); 
ylabel('Ez, V/m'); 
axis square xy; 

figure('Name','Ez_1D','NumberTitle','off'); 
plot(R1, Ez(p/2,:).*conj(Ez(p/2,:)));
xlabel('x, \mum'); 
ylabel('Ez, V/m'); 
axis square xy; 

figure('Name','E1_1D','NumberTitle','off'); 
plot(R1, Ez(p/2,:).*conj(Ez(p/2,:))+Ey(p/2,:).*conj(Ey(p/2,:))+Ex(p/2,:).*conj(Ex(p/2,:)));

figure('Name','Hx_1D','NumberTitle','off'); 
plot(R1, Hx(p/2,:).*conj(Hx(p/2,:)));

figure('Name','Hy_1D','NumberTitle','off'); 
plot(R1, Hy(p/2,:).*conj(Hy(p/2,:)));

figure('Name','Hz_1D','NumberTitle','off'); 
plot(R1, Hz(p/2,:).*conj(Hz(p/2,:)));

figure('Name','H1_1D','NumberTitle','off'); 
plot(R1, Hz(p/2,:).*conj(Hz(p/2,:))+Hy(p/2,:).*conj(Hy(p/2,:))+Hx(p/2,:).*conj(Hx(p/2,:)));

figure('Name','Sx','NumberTitle','off');
imagesc(R1, R1, Sx);
xlabel('x, \mum');
ylabel('y, \mum');
axis square xy;
colorbar

figure('Name','Sy','NumberTitle','off');
imagesc(R1, R1, Sy);
xlabel('x, \mum');
ylabel('y, \mum')
axis square xy;
colorbar

figure('Name','Sz','NumberTitle','off');
imagesc(R1, R3, Sz);
xlabel('x, \mum');
ylabel('z, \mum');
axis square xy;
colorbar

figure('Name','Sx_1D','NumberTitle','off'); 
plot(R1, Sx(p/2,:));
xlabel('x, \mum'); 
ylabel('Sx, a.u.'); 
axis square xy; 

figure('Name','Sy_1D','NumberTitle','off'); 
plot(R1, Sy(p/2,:));

figure('Name','Sz_1D','NumberTitle','off'); 
plot(R1, Sz(p/2,:));
xlabel('x, \mum'); 
ylabel('Sz, a.u.'); 
axis square xy; 

% figure('Name','x','NumberTitle','off'); 
% plot(graph_z, X1);
% xlabel('z, \lambda');
% ylabel('x\_min, W/m^2');
% axis square xy;
% 
% figure('Name','y','NumberTitle','off'); 
% plot(graph_z, Y1);
% xlabel('z, \lambda');
% ylabel('x\_max, a.u.');
% axis square xy;

