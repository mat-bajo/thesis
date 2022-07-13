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
        hmax = 1;
	case 2
        K = 1;
        A = [[0 0 1 0];[0 0 0 1];[-10-K 10 0 0];[5 -15 0 -0.25]];
        Ad = [[0 0 0 0];[0 0 0 0];[1 0 0 0];[0 0 0 0]];
        Bd = [0;0;K;0]; Cd = [1 0 0 0];
        [m,n] = size(Cd);
        hlist = {[0.1 2],[2 3],[3 4]};
        hminexpect = [0 2.673];
        hmaxexpect = [1.424 3.940];
        hmax = 4;
    otherwise
        error('Wrong choice of initial condition')
end

cmap = flip(parula(13));%hot;
Mlist = 0:5;%12;


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
        hminres = 0;
    else
        pocket = 1; % Stable
    end
    for ind = 1:(2*length(hmaxexpect) - 1*(pocket==1))
        hbound = hlist{ind};
        while hbound(2)-hbound(1) > step
            h = (hbound(1) + hbound(2))/2;
            res = posn(A,Bd,Cd,h,M);
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
            hmaxres = hbound(1);
            plot([M M],[hminres hmaxres],'Color',col,'LineWidth',4); hold on;
        end
        if pocket == 1
            hminres = hbound(2);
        end
        pocket = -pocket; % Become unstable || stable
    end
    if pocket == -1
        plot([M M],[hminres hmax],'Color',col,'LineWidth',4); hold on;
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



function res = posn(A,B,C,h,n)
    digits(32);%300);
    Ad = B*C;
    [nx,nz] = size(B);
    %% Calculation Gamman
    Mvpa = vpa([kron(A',eye(nx)) kron(Ad',eye(nx)); -kron(eye(nx),Ad') -kron(eye(nx),A')]);
    Mi = eye(2*nx^2)/Mvpa;
    eM = expm(h*Mvpa);
    Nvpa = [kron(A',eye(nx))+kron(eye(nx),A') kron(Ad',eye(nx)); eye(nx^2) zeros(nx^2)] + [kron(eye(nx),Ad') zeros(nx^2); zeros(nx^2) -eye(nx^2)]*expm(h*Mvpa);
    Ni = eye(2*nx^2)/Nvpa;
    % Recursive relations for Gamma
    Gamma = cell(n,1);
    Gamma{1} = Mi*(eM-eye(2*nx^2));
    Gamma{2} = -2/h*Mi*Gamma{1} + Mi*(eM+eye(2*nx^2));
    for k = 2:n-1
        Gamma{k+1} = Gamma{k-1} - 2*(2*k-1)/h*Mi*Gamma{k};
    end
    % Recursive relations for Gammabar
    Gammabar = cell(n);
    Gammabar{1,1} = Mi*(Gamma{1}-h*eye(2*nx^2));
    Gammabar{1,2} = -Mi*Gamma{2};
    Gammabar{2,2} = Mi*((2/h*Mi-eye(2*nx^2))*Gamma{2}-h/3*eye(2*nx^2));
    for i = 1:n
        for j = 1:n
            if j >= max(3,i)
                Gammabar{i,j} = Gammabar{i,j-2} + 2*(2*j-3)/h*Mi*Gammabar{i,j-1};
                if j == i
                    Gammabar{i,j} = Gammabar{i,j} - h/(2*i-1)*Mi;
                end
                if j == i+2
                    Gammabar{i,j} = Gammabar{i,j} + h/(2*i-1)*Mi;
                end
            else
                if j < i
                    Gammabar{i,j} = (-1)^(i+j)*Gammabar{j,i};
                end
            end
        end
    end
    % Calcul Pn
    U0 = reshape([eye(nx^2) zeros(nx^2)]*Ni*[-reshape(eye(nx),[],1); zeros(nx^2,1)],nx,nx);
    In = 1/h*kron(diag(2*(0:n-1)+1),eye(nx));
    Qn = zeros(nx,n*nx);
    for k = 0:n-1
        Qk = reshape([eye(nx^2) zeros(nx^2)]*Gamma{k+1}*Ni*[-reshape(eye(nx),[],1); zeros(nx^2,1)],nx,nx);
        Qn(:,k*nx+1:(k+1)*nx) = Qk'*Ad;%1/h*
    end
    Tn = zeros(n*nx);
    for i = 0:n-1
        for j = 0:n-1
            Tij = reshape([eye(nx^2) zeros(nx^2)]*Gammabar{i+1,j+1}*Ni*[-reshape(eye(nx),[],1); zeros(nx^2,1)],nx,nx);
            Tn(i*nx+1:(i+1)*nx,j*nx+1:(j+1)*nx) = Ad'*(Tij+(-1)^(i+j)*Tij')*Ad;%1/(h^2)*
        end
    end
    Pn = [U0 Qn; Qn' Tn+inv(In)];
    if min(double(real(eig(Pn)))) > 0
        res = 1; % Stable
    else
        res = 0; %Unstable
    end
end
