%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% GC rxn simulation %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc
%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%
%%% Time units are in hours %%%
%%% Simulation Parameters %%%
dt = 5e-3; % timestep
Dend = 21; % simulation termination day
N_T = 200; % number of T cells
NB0 = 250; % number of founder B cells
%%% T cell parameters %%%
% Variaility function form: P(A) = floor((max div + 1)((1-exp(-A/A0))
A0 = 5; % T cell selection parameter
%%% Fixed parameters %%%
tau_d = 8; % CC lifetime
r_ab = 3/4; % base antigen attempt rate
r_tb = 0.4;  % T-B cell encounter rate
tau_LD = 2.04; % LZ to DZ transition time
tau_DL = 5.88; % DZ of LZ transition time
tau_div = 5; % B cell division time
% mutation probabilities
pmut = 0.5; % probability of mutation upon division
pd = 0.3; % lethal mutation probability
pneg = 0.19; % negative mutation probability
ppos = 0.01; % positive mutation probability
p0 = 1-pd-pneg-ppos; % silent mutation probability
DelA = 1; % change in affinity due to mutation
% print/save options
Pflag = 1; % == 1 to print figures during simulation
Qpet = round(24/dt); % figure refresh duration
Rep = 1; % number of repetitions to run
%%%%%%%%%%%%%%%%%%%%% Execution %%%%%%%%%%%%%%%%%%%%%%
M = round(24*Dend/dt);
t = dt*(0:M)'; % time vector
ed_AN = -3.5:10.5; x_AN = ed_AN(1:end-1) + 0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mode: 0: fixed mutation rate, 1: variable mutation rate
mode_str = {'cnt','var'}; % mode
for mode = 0:1
    fprintf('mode %.0f... \n',mode);
    for runNUM = 1:Rep
        fprintf('rep %.0f out of %.0f \n',runNUM,Rep);
        % preallocate
        A_LZ = cell(Dend + 1,1); % affinity vector
        pop_all = cell(Dend + 1,1); % population vector
        CCA_all = cell(Dend + 1,1); % antigen vector
        B_all = cell(Dend+1,1);
        m_run = 0; % mutation ID
        Sel_Num = zeros(6,1);
        Sel_Out = cell(6,1);
        for i = 1:6
            Sel_Out{i} = zeros(length(ed_AN)-1,1);
        end
        %%%%%
        B_cells = initialize_Bcells(NB0,tau_d); % initialize B cells
        [pop,LZa,CC_A] = extract(B_cells); % GC properties
        B_all{1} = B_cells;
        A_LZ{1} = LZa; pop_all{1} = pop; CCA_all{1} = CC_A; % initial properties
        cc = 0; % day counter
        a0 = 0;
        %%%%%
        peakA = zeros(M,2);
        for i = 1:M
            % divison/mutation
            [B_cells,m_run, ANout] = DZ_division(B_cells, tau_div, pd, p0, pneg, ppos,pmut, dt,tau_DL, DelA, mode,m_run);
            for j = 1:6
                [stemp,~] = histcounts(ANout{j},ed_AN);
                Sel_Out{j} = Sel_Out{j} + stemp';
            end
            % CC --> CB differentiation
            B_cells = LZDZ_differ(B_cells,dt,tau_div);
            % TFH - B cell interaction
            [B_cells, Sel_temp] = TB_interaction(B_cells,N_T,r_tb,A0,dt,tau_LD);
            Sel_Num = Sel_Num + Sel_temp;
            % LZ B cell antigen collection
            B_cells = antigen_collection(B_cells,dt,r_ab,a0);
            % CB --> CC differentation
            B_cells = DZLZ_differ(B_cells,dt,tau_d);
            % CC apoptosis
            B_cells = cell_death(B_cells,dt);
            [pop,LZa,CC_A] = extract(B_cells); % GC properties
            a0 = mean(LZa);
            Amax = max(LZa); Nmax = length(find(LZa == Amax));
            peakA(i,:) = [Amax, Nmax];
            if ~mod(i,round(24/dt))
                cc = cc + 1;
                [pop,LZa,CC_A] = extract(B_cells); % GC properties
                B_all{cc + 1} = B_cells; % B cells
                A_LZ{cc+1} = LZa; pop_all{cc+1} = pop; CCA_all{cc+1} = CC_A; % GC properties
                if Pflag == 1
                    figure(1); clf;
                    y1 = pop_aver(pop_all(1:cc+1));
                    subplot(221);
                    plot(0:cc, y1,'ko-','markerface','k');
                    xlim([0,Dend]); ylim([0,8000])
                    [y2a,y2m,y2med] = A_aver(A_LZ(1:cc+1));
                    subplot(222);  hold on
                    plot(0:cc,y2a,'ko-','markerface','k');
                    plot(0:cc,y2m,'ro-','markerface','r');
                    plot(0:cc,y2med,'bx-');
                    xlim([0,Dend]); ylim([-1,5])
                    subplot(223);
                    ed3 = -0.5:(A0/3):(5*A0+0.5); x3 = ed3(1:end-1) + 0.5*(ed3(2)-ed3(1));
                    [y3,~] = histcounts(CC_A,ed3);
                    bar(x3,y3,'k');
                    yyaxis right
                    xr = linspace(0,5*A0+1,1e2);
                    Y = floor((6+1)*(1-exp(-((xr-0)/A0))));
                    Y(Y<0) = 0;
                    plot(xr,Y,'r.-','linewidth',2);
                    xlim([ed3(1),ed3(end)])
                    subplot(224);
                    ed4 = -5.5:5.5; x4 = ed4(1:end-1) + 0.5*(ed4(2)-ed4(1));
                    [y4,~] = histcounts(LZa,ed4);
                    bar(x4,y4,'k'); xlim([ed4(1),ed4(end)])
                    drawnow
                end
            end
        end
        out = {A_LZ,pop_all,CCA_all,B_all,Sel_Out,Sel_Num, peakA};
    end
end
%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%
function y = pop_aver(yin)
y = zeros(length(yin),1);
for i = 1:length(yin)
    y(i) = sum(yin{i});
end
end
function [ya,ym,ymed] = A_aver(yin)
ya = zeros(length(yin),1);
ym = ya;
ymed = ya;
for i = 1:length(yin)
    ya(i) = mean(yin{i});
    ym(i) = max(yin{i});
    ymed(i) = median(yin{i});
end
end
function Q = indexer(B_cells,phen) % pulls out relevant phenotype
Q = zeros(length(B_cells(:,1)),1);
cc = 0;
for i = 1:length(B_cells(:,1))
    if B_cells{i,1}(1) == phen
        cc = cc + 1;
        Q(cc) = i;
    end
end
Q = Q(1:cc);
end
function B_cells = initialize_Bcells(NB0,tau_d) % initialize B cells
B_cells = cell(NB0,2);
for i = 1:NB0
    % [identity, affinity, antigen, lifetime, num_div, org_num_div,org_aff]
    % identity: (0, CC), (1,CB), (2, CC->CB), (3, CB->CC)
    vs = 0;
    B_cells{i,1} = [0, vs, abs(randn(1)/2), tau_d,0,0,vs];
    B_cells{i,2} = nan(100,1); % mutation history
    B_cells{i,2}(1) = i; % assume all seeding cells are unique
end
end
%%%%%%%%%%%%%%%%%%%
function [B_cells,m_run, ANout] = DZ_division(B_cells, tau_div, pd, p0, pneg, ppos,pmut, dt,tau_DL, DelA, mode,m_run) % DZ division
% countdown mutation clock for DZ B cells and make mutations
ANout = cell(6,1);
for i = 1:6
    ANout{i} = nan(500,1);
end
ANvec = zeros(6,1);
% index of DZ B cells
Q = indexer(B_cells,1);
cc = length(Q);
if cc > 0.5
    for i = 1:cc
        a = B_cells(Q(i),:);
        b = {[],[]};
        a{1}(4) = a{1}(4) - dt;
        if a{1}(4) <= 0 % mutation upon division
            a{1}(4) = tau_div;
            a{1}(5) = round(a{1}(5) - 1);
            b = a; % division a --> a + b
            if rand(1) <= 1
                [a,m_run] = mutation(a,pd,p0,pneg,ppos,DelA,mode,pmut,m_run);
            end
            if rand(1) <= 1
                [b,m_run] = mutation(b,pd,p0,pneg,ppos,DelA,mode,pmut,m_run);
            end
        end
        % if done dividing, differentiate back to CC
        if ~isempty(a{1}) && a{1}(5) < 0.5
            a{1}(1) = 3;
            a{1}(4) = tau_DL;

            ANvec(a{1}(6)) = ANvec(a{1}(6)) + 1;
            ANout{a{1}(6)}(ANvec(a{1}(6))) = a{1}(2) - a{1}(7);
        end
        if ~isempty(b{1}) && b{1}(5) < 0.5
            b{1}(1) = 3;
            b{1}(4) = tau_DL;

            ANvec(b{1}(6)) = ANvec(b{1}(6)) + 1;
            ANout{b{1}(6)}(ANvec(b{1}(6))) = b{1}(2) - b{1}(7);
        end
        B_cells(Q(i),:) = a;
        B_cells(end+1,:) = b;
    end
    for i = 1:6
        ANout{i} = ANout{i}(1:ANvec(i));
    end
    B_cells = B_cells(~cellfun('isempty',B_cells(:,1)),:);
end
end
%%% LZ --> DZ migration %%%
function B_cells = LZDZ_differ(B_cells,dt,tau_div)
% countdown differentiation buffer, then send B cells to DZ
% index of LZ --> DZ differentiating B cells
Q = indexer(B_cells,2);
cc = length(Q);
% count down B cell differentiation clock
if length(cc) > 0.5
    for i = 1:cc
        a = B_cells{Q(i),1};
        a(4) = a(4) - dt;
        if a(4) <= 0 % finished differentiating
            a(1) = 1;
            a(4) = tau_div;
        end
        B_cells{Q(i),1} = a;
    end
end
end
%%% TB interaction %%%
function [B_cells, SelNum] = TB_interaction(B_cells,N_T,r_tb,A0,dt,tau_LD)
SelNum = zeros(6,1);
Amin = 2;
% index of LZ B cells
Q = indexer(B_cells,0); % LZ B cells
cc = length(Q);
if cc > 0.5
    % select LZ B cells
    inds = Q(randi(cc,N_T,1)); % randomly select a B cell
    for i = 1:N_T % for each TFH cell ...
        if logic_prob(r_tb,dt) % if they happen to run into a B cell ...
            ind = inds(i);
            obj = B_cells{ind,1};
            Y = floor((6+1)*(1-exp(-((obj(3)-Amin)/A0)))); % number of divisions
            if Y > 0.5
                srt = obj(2); % original affinity before entering DZ
                B_cells{ind,1} = [2,srt,0,tau_LD,Y,Y,srt];
                SelNum(Y) = SelNum(Y) + 1;
            end
        end
    end
end
end
%%% mutation machinery %%%
function [a, m_run] = mutation(a,pd,p0,pneg,ppos,DelA,mode,pmut,m_run)
% assign mutation
p1 = 0.6; p2 = 0.2; % mutation probability at D=1 and D=6
switch mode
    case 3 % no mutation
        y = 0;
    case 2
        m = (p2-p1)/(6-1);
        y = m*(a{1}(6)-1) + p1;
    case 1 % linear decay in mutation probability
        m = (p2-p1)/(6-1);
        y = m*(a{1}(6)-1) + p1;
    case 0
        y = pmut;
end
if rand(1) < y
    m_run = m_run + 1;
    pp = cumsum([pd,p0,pneg,ppos]);
    pp = pp/pp(end);
    rnum = rand(1);
    mode = find(pp > rnum, 1);
    switch mode
        case 1
            % B cell dies
            a{1} = [];
        case 2
            % nothing happens
            ind = find(isnan(a{2}),1);
            a{2}(ind) = m_run; % record mutation
        case 3
            % negative mutation
            a{1}(2) = a{1}(2) - DelA;
            ind = find(isnan(a{2}),1);
            a{2}(ind) = m_run; % record mutation
        case 4
            % positive mutation
            a{1}(2) = a{1}(2) + DelA;
            ind = find(isnan(a{2}),1);
            a{2}(ind) = m_run; % record mutation
    end
end
end
%%% antigen collection %%%
function B_cells = antigen_collection(B_cells,dt,r_ab,a0)
% index of LZ B cells
Q = indexer(B_cells,0); % LZ B cells
cc = length(Q);
if cc > 0.5
    for i = 1:length(Q)
        a = B_cells(Q(i),:);
        if a{1}(2) < 15
            SOS = r_ab; % ignore FDC interference
            a{1}(3) = a{1}(3) + 2*(1/(1+exp(a0-a{1}(2))))*logic_prob(SOS,dt); % stochastic antigen collection
            B_cells{Q(i),1}(3) = a{1}(3); % update antigen amount
        end
    end
end
end
%%% DZ --> LZ migration %%%
function B_cells = DZLZ_differ(B_cells,dt,tau_d)
% countdown differentiation buffer, then send B cells to DZ
% index of DZ --> LZ differentiating B cells
Q = indexer(B_cells,3);
cc = length(Q);
% count down B cell differentiation clock
if length(cc) > 0.5
    for i = 1:cc
        a = B_cells{Q(i),1};
        a(4) = a(4) - dt;
        if a(4) <= 0
            a(1) = 0;
            a(4) = tau_d;
        end
        B_cells{Q(i),1} = a;
    end
end
end
%%% Cell death %%%
function B_cells = cell_death(B_cells,dt)
% B cell apoptosis
% index of LZ B cells
Q = indexer(B_cells,0); % LZ B cells
cc = length(Q);
if length(cc) > 0.5
    for i = 1:cc
        a = B_cells{Q(i),1}(4);
        a = a - dt;
        B_cells{Q(i),1}(4) = a;
        if a <= 0
            B_cells{Q(i),1} = [];
        end
    end
end
B_cells = B_cells(~cellfun('isempty',B_cells(:,1)),:);
end
%%%%%%%%%%%%%%%%%
function [pop,LZa,CC_A] = extract(B_cells) % extract population and affinities
pop = zeros(4,1);
for i = 0:3
    Q = indexer(B_cells,i);
    pop(i+1) = length(Q);
end
LZa = zeros(length(B_cells),1);
for i = 1:length(B_cells)
    LZa(i) = B_cells{i,1}(2); % Affinity
end
temp = indexer(B_cells,0); % LZ B cells
CC_A = zeros(length(temp),1);
for i = 1:length(temp)
    CC_A(i) = B_cells{temp(i),1}(3);
end
end
%%% poisson %%%
function s = logic_prob(r,dt)
% returns 1 with probability r*dt, 0 otherwise
s = double(rand(1) < r*dt);
end