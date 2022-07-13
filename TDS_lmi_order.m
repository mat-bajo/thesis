% Lines
Klist = 10.^(-1:0.001:1);%10.^(-2:0.001:2);
hlist = 0.001:0.0005:3;%0.001:0.0005:10;
Knum = length(Klist);
hnum = length(hlist);

Zu = NaN(Knum,hnum);
for Kind = 1:Knum
    K = Klist(Kind);
    A = [[0 0 1 0];[0 0 0 1];[-10-K 10 0 0];[5 -15 0 -0.25]];
    B = [0;0;K;0]; C = [1 0 0 0];
    [nx,nz] = size(B);

    M = allmargin(ss(A,-B,C,0));
    Del = M.DelayMargin;
    Del = [(M.PhaseMargin*pi/180+2*pi)./M.PMFrequency Del]; 
    Del = [(M.PhaseMargin*pi/180+4*pi)./M.PMFrequency Del]; 
    Del = [(M.PhaseMargin*pi/180+6*pi)./M.PMFrequency Del];
    Del = [(M.PhaseMargin*pi/180+8*pi)./M.PMFrequency Del];
    while not(isempty(Del))
        for hind = 1:hnum
            h = hlist(hind);
            if h > Del(end) && h < Del(end)+0.001
                Zu(Kind,hind) = 1;
            end
        end
        Del = Del(1:end-1);
    end
end

Sol = NaN(Knum,hnum);
for Kind = 1:Knum
    for hind = 1:hnum
        if Zu(Kind,hind) == 1
            Sol(Kind,hind) = 1;
        end
    end
end

[X,Y] = meshgrid(log10(Klist),hlist);
plot3(X,Y,Sol','.k','Linewidth',0.1);hold on; view(2); 
set(gca,'FontSize',20)

% Parameters
Klist = 10.^(-1:0.05:1);
hlist = 0:0.05:3;%[0.1 0.5 1 2 5 10];%
Knum = length(Klist);
hnum = length(hlist);
Z = zeros(Knum,hnum);

for Kind = 1:Knum
    % System
    K = Klist(Kind);
    A = [[0 0 1 0];[0 0 0 1];[-10-K 10 0 0];[5 -15 0 -0.25]];
    B = [0;0;K;0]; C = [1 0 0 0];
    for hind = 1:hnum
       h = hlist(hind);
       if log10(nast_lmi(A,B,C,h)) < 50
           Z(Kind,hind) = log10(nast_lmi(A,B,C,h));
       else
           Z(Kind,hind) = 50;
       end
    end
    disp(['Step:' int2str(Kind) '/' int2str(Knum)]);
end

[X,Y] = meshgrid(log10(Klist),hlist);
surf(X,Y,Z','EdgeColor','none'); view(2); %,'FaceColor','k','FaceAlpha',.3,'EdgeAlpha',.3)
alpha 0.75
xlabel('Parameter $log_{10}(\lambda)$','Interpreter','Latex')
ylabel('Delay $h$','Interpreter','Latex')
zlabel('Order $N^{\ast}$','Interpreter','Latex')
hcb = colorbar('TickLabels',{'10^0','10^{10}','10^{20}','10^{30}','10^{40}','10^{50}'},'Ticks', [1 10:10:50]);
%hcb.Title.String = 'Order N^\ast';
ylabel(hcb,'Order $N^{\ast}$','Interpreter','Latex')
grid on; set(gca, 'fontsize', 22);

function nast = nast_lmi(A,B,C,h)
% Calcultaiton of order nast
    
    [nx,nz] = size(B);
    
    % Lyapunov matrix function f 
    Ad = B*C;
    
    %digitsOld = digits(32);
    Mcal = -[kron(A',eye(nx)) kron(Ad',eye(nx)); -kron(eye(nx),Ad') -kron(eye(nx),A')];%vpa()
    Ncal = [kron(A',eye(nx))+kron(eye(nx),A') kron(Ad',eye(nx)); eye(nx^2) zeros(nx^2)] + [kron(eye(nx),Ad') zeros(nx^2); zeros(nx^2) -eye(nx^2)]*expm(-h*Mcal);
    Ni = eye(2*nx^2)/Ncal;
    
    % Frobenius norm of several matrices
    M1 = Mcal^2*Ni;
    M2 = Mcal^4*Ni;
    Mnorm = sqrt(max(eig(Mcal'*Mcal)));
    M1norm = sqrt(max(eig(M1'*M1)));
    M2norm = sqrt(max(eig(M2'*M2)));
    Anorm = sqrt(max(eig(A'*A)));
    Cnorm = sqrt(max(eig(C'*C)));
    Bnorm = sqrt(max(eig(B'*B)))*Cnorm;

    % Decay rates
    rho1 = sqrt(nx)*(pi/2)^(3/2)*exp(h*Mnorm)*M1norm*h^2*Bnorm;
    rho2 = 1/2*sqrt(nx)*(pi/2)^(3/2)*exp(h*Mnorm)*M2norm*h^3*Bnorm;
    rho3 = sqrt(2*pi)*(1+sqrt(nx)*pi*h*exp(h*Mnorm)*M1norm)*h*Bnorm;
    
    % Calcul de l'ordre minimal requis nast
    nast = 5+ceil((1+h^2)^2*Bnorm^2*(rho1*(4/(1+h^2)+Anorm+Bnorm)+rho2+2*rho3)^2);
    %digits(digitsOld);
    
end

