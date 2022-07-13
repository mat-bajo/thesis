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
Klist = 10.^(-1:0.05:2);
hlist = 0.01:0.01:10;

for Kind = 1:length(Klist)
    % System
    K = Klist(Kind);
    A = [[0 0 1 0];[0 0 0 1];[-10-K 10 0 0];[5 -15 0 -0.25]];
    B = [0;0;K;0]; C = [1 0 0 0];
    Ad = B*C;
    [nx, nz]=size(B);
    for hind = 1:length(hlist)
       h = hlist(hind);
       Mcal = h*[kron(A',eye(nx)) kron(Ad',eye(nx)); -kron(eye(nx),Ad') -kron(eye(nx),A')];
       Ncal = [kron(A',eye(nx))+kron(eye(nx),A') kron(Ad',eye(nx)); eye(nx^2) zeros(nx^2)] + [kron(eye(nx),Ad') zeros(nx^2); zeros(nx^2) -eye(nx^2)]*expm(h*Mcal);
       Ncalinv = inv(Ncal);
        rho0 = exp(normM(Mcal))*norm(Ncalinv);
        kappa1 = h*normM(Ad)*rho0;
        kappa2 = h^2*normM(Ad)*rho0;
        r = normM(A)+normM(Ad);
        mu = h*r/2;
        rho0 = 1;
        a = kappa2;
        b = kappa1+kappa2;
        c = 1/(2*r);
        epsilon = Eps(a,b,c);
        order = Nasteps(epsilon,mu,rho0);
        if order > 300 || rho0 == Inf || a == Inf 
            Z(Kind,hind) = 300;
        else
            Z(Kind,hind) = order;
        end
    end
    disp(['Step:' int2str(Kind) '/' int2str(length(Klist))]);
end

[X,Y] = meshgrid(log10(Klist),hlist);
surf(X,Y,Z','EdgeColor','none'); view(2); hold on; %,'FaceColor','k','FaceAlpha',.3,'EdgeAlpha',.3)
alpha 0.75
ylabel('Delay $h$','Interpreter','Latex')
colormap(flip(parula))
hcb = colorbar('TickLabels',{'50','100','150','200','250','>300'},'Ticks',[50 100 150 200 250 300]);
xlabel('Parameter $log_{10}(k)$','Interpreter','Latex')
ylabel(hcb,'Order $n$','Interpreter','Latex')%such that LMIs hold
grid on; set(gca, 'fontsize', 22);

function norm = normM(M)
    norm = sqrt(max(eig(M'*M)));%trace
end

function out = Eps(a,b,c)
    % Solve 1-bx-cx^2<0
    out = c / (b+sqrt(b^2+a*c));
    %-b/a + sqrt(b^2+a*c)/a;
end

function nast = Nasteps(eps,r,rho0)
% Theorem 1: Calculation of N1(eps)
    rc = floor(r);
    mu = exp(rc)/(2*r^2*rho0);
    nast = double(abs(ceil(2 + r*exp(1+lambertw(-log(mu*eps)/(r*exp(1)))))));
end

