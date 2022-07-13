clear all; close all;

Mmax = 12;
num = input('Choose the system (1 or 2):');
disp('Creation of the system...')

switch num
    case 1
        A = 1;
        Ad = -2;
        Bd = -2; Cd = 1;
        [m,n] = size(Cd);
        h = 0.3;
        Rayon0 = 5;
        Rayon = 0.01;
        Mlist = 0:2:4;%1:3;
        xb = [-1.9-0.05,-1.9+0.05]; yb = [-2.05-0.05,-2.05+0.05];
        posl = [0.72 0.28 0.11 0.22];
        numlist = 2;
    case 2
        K=1;
        A = [[0 0 1 0];[0 0 0 1];[-10-K 10 0 0];[5 -15 0 -0.25]];
        Ad = [[0 0 0 0];[0 0 0 0];[1 0 0 0];[0 0 0 0]];
        Bd = [0;0;K;0]; Cd = [1 0 0 0];
        [m,n] = size(Cd);
        h = 0.3;
        Rayon0 = 5;
        Rayon = 0.01;
        Mlist = 0:3;
        xb = [-0.13-0.02,-0.13+0.02]; yb = [-2.25-0.05,-2.25+0.05];
        posl = [0.74 0.28 0.11 0.22];
        numlist = [2 4];
    otherwise
        error('Wrong choice of initial condition')
end


E0 = cell(Mmax+1,1);
E1 = cell(Mmax+1,1);

for M = 0:Mmax
    %Useful matrices
    I0M = (-1).^(0:M-1)';
    I1M = ones(M,1);
    U0M = I0M*I0M';
    U1M = I1M*I1M';
    DM = diag((2*(0:M-1)+1))/h;
    LM = tril(I1M*I1M'-I0M*I0M',-1);
    % Delay
    AM = -(LM+U0M)*DM;
    B1M = I1M;
    B0M = I0M;
    C1M = I1M'*DM;
    C0M = I0M'*DM;
    % System (x,XN)
    At = [A Bd*kron(C0M,eye(m)); kron(B1M,eye(m))*Cd kron(AM,eye(m))];
    lambda = eig(At);                       %approximated eigenvalues
    [~,ind] = sort(abs(lambda));%,'descend'); %ordered by decreasing real part
    lambda = lambda(ind);
    E0{M+1} = [real(lambda) imag(lambda)];
    
    %Useful matrices
    I0M = (-1).^(0:M-1)';
    I1M = ones(M,1);
    U0M = I0M*I0M';
    U1M = I1M*I1M';
    DM = diag((2*(0:M-1)+1))/h;
    LM = tril(I1M*I1M'-I0M*I0M',-1);
    % Delay
    AM = (LM'-U1M-(-1)^(M-1)*I0M*I1M')*DM;
    BM = I1M + (-1)^(M-1)*I0M;
    CM = (I0M' + (-1)^(M-1)*I1M')*DM;
    DM = (-1)^M;
    B0M = I0M;
    % System (x,XN)
    At = [A+Bd*kron(DM,eye(m))*Cd Bd*kron(CM,eye(m)); kron(BM,eye(m))*Cd kron(AM,eye(m))];
    lambda = eig(At);   %approximated eigenvalues
    [~,ind] = sort(abs(lambda));%,'descend'); %ordered by decreasing real part
    lambda = lambda(ind);
    E1{M+1} = [real(lambda) imag(lambda)];
end

M = 50;
%Useful matrices
I0M = (-1).^(0:M-1)';
I1M = ones(M,1);
U0M = I0M*I0M';
U1M = I1M*I1M';
DM = diag((2*(0:M-1)+1))/h;
LM = tril(I1M*I1M'-I0M*I0M',-1);
% Delay
AM = -(LM+U0M+(-1)^(M-1)*I0M*I1M')*DM;%(LM'-U1M-(-1)^(M-1)*I0M*I1M')*DM;% 
BM = I1M + (-1)^(M-1)*I0M;
CM = (I0M' + (-1)^(M-1)*I1M')*DM;
DM = (-1)^M;
B0M = I0M;
% System (x,XN)
At = [A+Bd*kron(DM,eye(m))*Cd Bd*kron(CM,eye(m)); kron(BM,eye(m))*Cd kron(AM,eye(m))];
lambda = eig(At);   %approximated eigenvalues
[~,ind] = sort(abs(lambda));%,'descend'); %ordered by decreasing real part
lambda = lambda(ind);
Eex = [real(lambda(1:n+m*Mmax)) imag(lambda(1:n+m*Mmax))];

% Display
cmap = flip(parula(13));%hot;
figure(1)
XCentre = 0;
YCentre = 0;
VTheta = 0:0.1:180;
XCercle = XCentre + Rayon0 * cos(VTheta);
YCercle = YCentre + Rayon0 * sin(VTheta);
plot(XCercle, YCercle, 'k-','DisplayName','$\mathcal{B}(0,R)$','LineWidth',0.5); hold on;
%E
for M=0:2:Mmax
    col = cmap(M+1,:);
    %scatter(E0{M+1}(:,1),E0{M+1}(:,2),200,col,'+','LineWidth',1.5,'DisplayName',['N=' int2str(M)]);hold on;
    %scatter(E1{M+1}(:,1),E1{M+1}(:,2),200,col,'x','LineWidth',1.5,'DisplayName',['N=' int2str(M)]);hold on;
    for indice=1:n+m*M
        plot(E0{M+1}(indice,1),E0{M+1}(indice,2),'+','Color',col,'MarkerSize',15,'LineWidth',2); hold on;
        plot(E1{M+1}(indice,1),E1{M+1}(indice,2),'x','Color',col,'MarkerSize',15,'LineWidth',2); hold on;
    end
end
scatter(Eex(1:11,1),Eex(1:11,2),50,'ko','filled','DisplayName','Expected')%16
xlabel('Real part','Interpreter','Latex'); ylabel('Imaginary part','Interpreter','Latex');
grid on; grid minor; 
set(gca, 'fontsize', 22);
xlim([-120 120]);ylim([-100 100]);
set(gcf,'Position',[100 100 1000 420]);
colormap(cmap);
hcb1 = colorbar('Direction','reverse','Ticks',((0:2:12)+0.5)/13,'TickLabels',{'0','2','4','6','8','10','12'});
ylabel(hcb1,'Order $n$','Interpreter','Latex')
p = gca;
% Calculate x,y points of zoomPlot
pos = [0.52 0.3 0.3 0.6]; %[left, bottom, width, height] specifies the location and size of the side of the zoom box
x1 = (pos(1)-p.Position(1))/p.Position(3)*diff(p.XLim)+p.XLim(1);
x2 = (pos(1)+pos(3)-p.Position(1))/p.Position(3)*diff(p.XLim)+(p.XLim(1));
y1 = (pos(2)-p.Position(2))/p.Position(4)*diff(p.YLim)+p.YLim(1);
y2 = ((pos(2)+pos(4)-p.Position(2))/p.Position(4))*diff(p.YLim)+p.YLim(1);
% Plot lines connecting zoomPlot to original plot points
plot([0 x1], [Rayon0 y2], ':k'); % Line to vertex 1
%plot([Rayon0 x2], [Rayon0 y2], ':k'); % Line to vertex 2
%plot([Rayon0 x2], [-Rayon0 y1], ':k'); % Line to vertex 4
plot([0 x1], [-Rayon0 y1], ':k'); % Line to vertex 3
% Plot zoomPlot and change axis
z = axes('position',pos);
%box on 
plot(XCercle, YCercle, 'k-','DisplayName','$\mathcal{B}(0,R)$'); hold on;
for num1 = numlist
    XCentre = Eex(num1,1);
    YCentre = Eex(num1,2);
    VTheta = 0:0.1:180;
    XCercle = XCentre + Rayon * cos(VTheta);
    YCercle = YCentre + Rayon * sin(VTheta);
    plot(XCercle, YCercle, 'k-','DisplayName','$\mathcal{B}(s^{\ast},r)$'); hold on;
end
%E
for M=Mlist
    col = cmap(M+1,:);%cmap(floor(64*(M-min(Mlist)+1)/(max(Mlist)-min(Mlist)+1)),:);
    %scatter(E0{M+1}(:,1),E0{M+1}(:,2),200,col,'+','LineWidth',1.5,'DisplayName',['N=' int2str(M)]);%hold on;
    %scatter(E1{M+1}(:,1),E1{M+1}(:,2),200,col,'x','LineWidth',1.5,'DisplayName',['N=' int2str(M)])
    for indice=1:n+m*M
        plot(E0{M+1}(indice,1),E0{M+1}(indice,2),'+','Color',col,'MarkerSize',15,'LineWidth',2); hold on;
        plot(E1{M+1}(indice,1),E1{M+1}(indice,2),'x','Color',col,'MarkerSize',15,'LineWidth',2); hold on;
    end
end
scatter(Eex(1:10,1),Eex(1:10,2),50,'ko','filled','DisplayName','Expected')
grid on; xlim([-Rayon0-1,Rayon0+1]); ylim([-Rayon0-1,Rayon0+1]); 
xticks([-Rayon0 Rayon0]);%xticklabels({'-5','5'})
yticks([-Rayon0 Rayon0]);%yticklabels({'-5','5'})
%t= text(-5*Rayon0/6,Rayon0/6,'B(0,R=5)');
%set(t,'Interpreter','tex', 'fontsize', 18)
set(gca, 'fontsize', 18);
%set(gca,'box','off');
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
% Get current axis position and limits
p = gca;
% Calculate x,y points of zoomPlot
pos = posl; %[left, bottom, width, height] specifies the location and size of the side of the zoom box
x1 = (pos(1)-p.Position(1))/p.Position(3)*diff(p.XLim)+p.XLim(1);
x2 = (pos(1)+pos(3)-p.Position(1))/p.Position(3)*diff(p.XLim)+(p.XLim(1));
y1 = (pos(2)-p.Position(2))/p.Position(4)*diff(p.YLim)+p.YLim(1);
y2 = ((pos(2)+pos(4)-p.Position(2))/p.Position(4))*diff(p.YLim)+p.YLim(1);
% Plot lines connecting zoomPlot to original plot points
plot([xb(1) x1], [yb(1) y2], ':k'); % Line to vertex 1
%plot([xb(2) x2], [yb(2) y2], ':k'); % Line to vertex 2
%plot([xb(2) x2], [yb(2) y1], ':k'); % Line to vertex 4
plot([xb(1) x1], [yb(1) y1], ':k'); % Line to vertex 3
% Plot zoomPlot and change axis
z = axes('position',pos);
box on 
for num1 = numlist
    XCentre = Eex(num1,1);
    YCentre = Eex(num1,2);
    VTheta = 0:0.1:180;
    XCercle = XCentre + Rayon * cos(VTheta);
    YCercle = YCentre + Rayon * sin(VTheta);
    plot(XCercle, YCercle, 'k-','DisplayName','$\mathcal{B}(s^{\ast},r)$'); hold on;
end
for M=Mlist
    col = cmap(M+1,:);%cmap(floor(64*(M-min(Mlist)+1)/(max(Mlist)-min(Mlist)+1)),:);
    %scatter(E0{M+1}(:,1),E0{M+1}(:,2),200,col,'+','LineWidth',1.5,'DisplayName',['N=' int2str(M)]);hold on;
    %scatter(E1{M+1}(:,1),E1{M+1}(:,2),200,col,'x','LineWidth',1.5,'DisplayName',['M=' int2str(M)])
    for indice=1:n+m*M
        plot(E0{M+1}(indice,1),E0{M+1}(indice,2),'+','Color',col,'MarkerSize',15,'LineWidth',2); hold on;
        plot(E1{M+1}(indice,1),E1{M+1}(indice,2),'x','Color',col,'MarkerSize',15,'LineWidth',2); hold on;
    end
end
scatter(Eex(1:10,1),Eex(1:10,2),50,'ko','filled','DisplayName','Expected')
grid on; grid minor; xlim(xb); ylim(yb)
set(gca, 'fontsize', 15);
pause(0.1)

figure(2)
for M=0:2:Mmax
    col = cmap(M+1,:);
    R = Eex(1:n+m*M,1).^2 + Eex(1:n+m*M,2).^2;
    r0 = (E0{M+1}(:,1) - Eex(1:n+m*M,1)).^2+(E0{M+1}(:,2) - Eex(1:n+m*M,2)).^2;
    r1 = (E1{M+1}(:,1) - Eex(1:n+m*M,1)).^2+(E1{M+1}(:,2) - Eex(1:n+m*M,2)).^2;
    %scatter(10*log10(R(indice)),10*log10(r0(indice)),200,col,'+','Linewidth',2); hold on;
    %scatter(10*log10(R(indice)),10*log10(r1(indice)),200,col,'x','Linewidth',2); hold on;
    for indice=1:n+m*M
        loglog(sqrt(R(indice)),sqrt(r0(indice)),'+','Color',col,'MarkerSize',15,'LineWidth',2); hold on;
        loglog(sqrt(R(indice)),sqrt(r1(indice)),'x','Color',col,'MarkerSize',15,'LineWidth',2); hold on;
    end
end
xlabel('$|s^{\ast}|$','Interpreter','Latex'); ylabel('$|{s^n}-s^{\ast}|$','Interpreter','Latex');
xticks([1 10 100]); xticklabels({'10^{0}','10^{1}','10^{2}'}); xlim([1 100])
yticks([10^(-15) 10^(-10) 10^(-5) 1 10^5]); yticklabels({'10^{-15}','10^{-10}','10^{-5}','10^{0}','10^5'}); ylim([10^(-15) 10^5])
grid on; set(gca, 'fontsize', 22); %set(gca,'TickLabelInterpreter', 'Latex');

colormap(cmap);
hcb2 = colorbar('Direction','reverse','Ticks',((0:2:12)+0.5)/13,'TickLabels',{'0','2','4','6','8','10','12'});
ylabel(hcb2,'Order $n$','Interpreter','Latex')