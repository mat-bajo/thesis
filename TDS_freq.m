clear all;
close all;

num = input('Choose the system (1 or 2):');
disp('Creation of the system...')
nbDecimal = 2; % Precision of the path-following algorithm
step = 10^(-nbDecimal); % Pas

Emax1 = [1 1.2596 1.3840 1.4690 1.5345 1.5883 1.6342 1.6744 1.7102 1.7426 1.7722 1.7995 1.8249 1.8485 1.8708 1.8918 1.9116 1.9305 1.9485 1.9657 1.9822 1.9980];
Emax2 = [2 2.1263 2.2217 2.2987 2.3636 2.4200 2.4700 2.5151 2.5562 2.5940 2.6290 2.6618 2.6925 2.7215 2.7489 2.7749 2.7998 2.8235 2.8462 2.8680 2.8889];
%floor(1./Emax1*10000)/10000
%floor(1./Emax2*10000)/10000
switch num
    case 1
        A = 1;
        Ad = -2;
        Bd = -2; Cd = 1;
        [m,n] = size(Cd);
        hlist = step:step:3;
        hminexpect = [];
        hmaxexpect = 0.604;
	case 2
        A = [[0 0 1 0];[0 0 0 1];[-4 0 0 0];[0 -16 0 0]];
        Ad = [[0 0 0 0];[0 0 0 0];[0 0 0 -1];[0 0 1 0]];
        Bd = [0 0;0 0;0 -1;1 0]; Cd = [0 0 1 0; 0 0 0 1];
        [m,n] = size(Cd);
        hlist = step:step:5;
        hminexpect = [0.4108 2.054];% 3.697];
        hmaxexpect = [0.7509 2.252];% 3.754];
    otherwise
        error('Wrong choice of initial condition')
end

cmap = flip(parula(13));%hot;
Mlist = 1:12;

%% Pade (N-1,N)

%% Search of the stability on hlist
disp('===================================');
disp('Direct stability test on AN');
disp('===================================');

%% Search of min-max delay h allowable
for M = Mlist
    col = cmap(M+1,:);%cmap(floor(64*(M+1)/(Mlist(end)+1)),:);
    disp('===================================');
    disp(strcat(['Is Pade(' int2str(M-1) ',' int2str(M) ') A' int2str(M) ' stable']));
    disp('===================================');
    if max(real(eig(A+Ad))) >= 0
        pocket = -1; % Unstable
        res = 0;
    else
        pocket = 1; % Stable
        res = 1;
        hminres = step;
    end
    hmin = 0;
    %Useful matrices
    I0M = (-1).^(0:M-1)';
    I1M = ones(M,1);
    U0M = I0M*I0M';
    U1M = I1M*I1M';
    LM = tril(I1M*I1M'-I0M*I0M',-1);
    B1M = I1M;
    for h = hlist
        % Matrix delay dependent
        DeltaM = diag((2*(0:M-1)+1))/h;
        % Delay
        AM = -(LM+U0M')*DeltaM;
        C0M = I0M'*DeltaM;
        At = [A Bd*kron(C0M,eye(m)); kron(B1M,eye(m))*Cd kron(AM,eye(m))];
        Bt = [Bd; -kron(I0M,eye(m))];
        Ct = [Cd zeros(m,m*M)];
        sys = ss(At,Bt,Ct,zeros(m));
        norm = hinfnorm(sys);
        if res == 0 && norm < 1/Emax1(M+1)
            res = 1;
            hminres = h;
        end
        if res == 1 && norm >= 1/Emax1(M+1)
            res = 0;
            hmaxres = h-step;
        end
        if (pocket == 1 && max(real(eig(At))) >= 0) || (pocket == -1 && max(real(eig(At))) < 0)
            hmax = h-step;
            if num == 1
                if pocket == 1 && hmaxres ~= 0
                    plot([M-0.1 M-0.1],[0 hmaxres],'Color',col,'Marker','+','MarkerIndices',2,'MarkerSize',15,'LineWidth',4); hold on;
                end
            else
                if pocket == 1 && hminres ~= hmaxres
                    plot([M-0.1 M-0.1],[hminres hmaxres],'Color',col,'Marker','+','MarkerSize',15,'LineWidth',4); hold on;
                end
            end
            pocket = -pocket; % Become unstable || stable
            hmin = h;
        end
        hmax = h;
    end
    disp('===================================');
    disp(strcat(['Is Pade(' int2str(M) ',' int2str(M) ') A' int2str(M) ' stable']));
    disp('===================================');
    if max(real(eig(A+Ad))) >= 0
        pocket = -1; % Unstable
        res = 0;
        hminres = 0;
        hmaxres = 0;
    else
        pocket = 1; % Stable
        res = 1;
        hminres = step;
    end
    hmin = 0;
    hmax = 0;
    % Delay
    BM = I1M + (-1)^(M-1)*I0M;
    DM = (-1)^M;
    for h = hlist
        % Matrix delay dependent
        DeltaM = diag((2*(0:M-1)+1))/h;
        % Delay
        AM = (LM'-U1M-(-1)^(M-1)*I0M*I1M')*DeltaM;%LM
        CM = (I0M' + (-1)^(M-1)*I1M')*DeltaM;
        At = [A+Bd*kron(DM,eye(m))*Cd Bd*kron(CM,eye(m)); kron(BM,eye(m))*Cd kron(AM,eye(m))];
        Bt = [Bd; -kron(I0M,eye(m))];
        Ct = [Cd zeros(m,m*M)];
        sys = ss(At,Bt,Ct,zeros(m));
        norm = hinfnorm(sys);
        if res == 0 && norm < 1/Emax2(M+1)
            res = 1;
            hminres = h;
        end
        if res == 1 && norm >= 1/Emax2(M+1)
            res = 0;
            hmaxres = h-step;
        end
        if (pocket == 1 && max(real(eig(At))) >= 0) || (pocket == -1 && max(real(eig(At))) < 0)
            hmax = h-step;
            if num ==1
                if pocket == 1 && hmaxres ~= 0
                    plot([M+0.1 M+0.1],[0 hmaxres],'Color',col,'Marker','x','MarkerIndices',2,'MarkerSize',15,'LineWidth',4); hold on;
                end
            else
                if pocket == 1 && hminres ~= hmaxres
                    plot([M+0.1 M+0.1],[hminres hmaxres],'Color',col,'Marker','x','MarkerSize',15,'LineWidth',4); hold on;
                end
            end
            pocket = -pocket; % Become unstable || stable
            hmin = h;
        end
        hmax = h;
    end
end

for h = [hminexpect hmaxexpect]
    stop = Mlist(end)+0.1;
    plot([-0.1:0.1:stop],h*ones(size(-0.1:0.1:stop)),'k.','Linewidth',2); hold on;
end
xlabel('Order $n$','Interpreter','Latex'); ylabel('Delay $h$','Interpreter','Latex');
xlim([-0.1 Mlist(end)+0.1]); ylim([0 ceil(hmaxexpect(end))])
grid on; set(gca, 'fontsize', 22);
%colormap(cmap);
%colorbar('Direction','reverse','Ticks',((0:2:8)+1)/11,'TickLabels',{'N=0','N=2','N=4','N=6','N=8','N=10'});