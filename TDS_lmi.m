clear all;
close all;

num = input('Choose the system (1 or 2):');
disp('Creation of the system...')
nbDecimal = 2; % Precision of the path-following algorithm
step = 10^(-nbDecimal); % Pas

switch num
    case 1
        A = 1;
        Ad = -2;
        Bd = -2; Cd = 1;
        [m,n] = size(Cd);
        hlist = {[0.5 1]};
        hminexpect = [];
        hmaxexpect = 0.604;
	case 2
        K = 1;
        A = [[0 0 1 0];[0 0 0 1];[-10-K 10 0 0];[5 -15 0 -0.25]];
        Ad = [[0 0 0 0];[0 0 0 0];[1 0 0 0];[0 0 0 0]];
        Bd = [0;0;K;0]; Cd = [1 0 0 0];
        [m,n] = size(Cd);
        hlist = {[0.1 2],[2 3],[3 4]};
        hminexpect = [0 2.673];
        hmaxexpect = [1.424 3.940];
    otherwise
        error('Wrong choice of initial condition')
end

cmap = flip(parula(13));%hot;
Mlist = 1:12;


%% Search of the stability on hlist
disp('===================================');
disp('LMI stability test');
disp('===================================');

%% Search of min-max delay h allowable
for M = Mlist
    disp(['Order N=' int2str(M)]);
    col = cmap(M+1,:);
    if max(real(eig(A+Ad))) >= 0
        pocket = -1; % Unstable
    else
        pocket = 1; % Stable
        hminres = 0;
    end
    for ind = 1:(2*length(hmaxexpect) - 1*(pocket==1))
        hbound = hlist{ind};
        while hbound(2)-hbound(1) > step
            h = (hbound(1) + hbound(2))/2;
            res = lmin(A,Bd,Cd,h,M);
            if pocket == -1
                if res == 0
                    hbound(1) = h;
                else
                    hbound(2) = h;
                end
            end
            if pocket == 1
                if res == 0
                    hbound(2) = h;
                else
                    hbound(1) = h;
                end
            end
        end
        if pocket == -1
            hminres = hbound(2);
        end
        if pocket == 1
            hmaxres = hbound(1);
            plot([M M],[hminres hmaxres],'Color',col,'LineWidth',4); hold on;
        end
        pocket = -pocket; % Become unstable || stable
    end
end

for h = [hminexpect hmaxexpect]
    stop = Mlist(end)+0.1;
    plot(-0.1:0.1:stop,h*ones(size(-0.1:0.1:stop)),'k.','Linewidth',2); hold on;
end
xlabel('Order $n$','Interpreter','Latex'); ylabel('Delay $h$','Interpreter','Latex');
xlim([-0.1 Mlist(end)+0.1]); ylim([0 ceil(hmaxexpect(end))])
grid on; set(gca, 'fontsize', 22);
%colormap(cmap);
%colorbar('Direction','reverse','Ticks',((0:2:8)+1)/11,'TickLabels',{'N=0','N=2','N=4','N=6','N=8','N=10'});



function res = lmin(A,B,C,h,n)

    [nx,nz] = size(B);
    nlist = (0:n-1)';
    Ian = kron((-1).^nlist,eye(nz)); % Boundary th=0
    Ibn = kron(1.^nlist,eye(nz)); % Boundary th=1
    In = kron(diag(2*nlist+1),eye(nz));
    Ln = tril(Ibn*Ibn'-Ian*Ian'); % Derivation

    %% Model A
    An = [A B*Ian'; 1/h*In*Ibn*C 1/h*In*(-Ln-Ian*Ian')];
    Bn = [B; -1/h*In*Ian];
    Cn = [C zeros(nz,nz*n)];
    Can = [zeros(nz,nx) -Ian'];

    % Lyapunov analysis
    %% LMI in S,R
    % Yalmip Solver
    Pn = sdpvar(nx+n*nz,nx+n*nz);
    S = sdpvar(nz,nz);
    R = sdpvar(nz,nz);
    %Construction of an LMI
    % Positivity of V
    Phi = Pn + [zeros(nx,nx) zeros(nx,n*nz); zeros(n*nz,nx) fIn(h*S,n)+fIn(h*R/2,n)+fIn(eye(nz),n)*fJn(h*R/2,n)*fIn(eye(nz),n)];
    % Negativity of \dot{V}
    Psi0 = [Pn*An+An'*Pn+Cn'*(S+R)*Cn-Can'*S*Can Pn*Bn+Can'*S; Bn'*Pn+S*Can -S];
    Psi = Psi0 - [zeros(nx,nx+(n+1)*nz); zeros(n*nz,nx) fIn(R,n) zeros(n*nz,nz); zeros(nz,nx+(n+1)*nz)];
    constraints = [Phi >= 1e-5, S >= 1e-5, R >= 1e-5, Psi <= -1e-5];
    options = sdpsettings('verbose',0,'solver','sdpt3');
    optimize(constraints,[],options);
    pres = checkset(constraints);
    if sum(pres > 0) < length(pres)
        res = 0;
    else
        res = 1;
    end
end

function Mn = fIn(M,n)
    nlist = (0:n-1)';
    Dn = diag(1./(2*nlist+1));
    Mn = kron(Dn,M);
end

function Mn = fJn(M,n)
    nlist = (1:n-1)';
    Dmn = diag(nlist,-1);
    Dpn = diag(nlist,1);
    Mn = kron(Dmn,M)+kron(Dpn,M);
end