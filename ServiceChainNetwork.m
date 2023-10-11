% Flow model simulator for the joint routing/compuation with general convex cost.
% Created by Jinkun, Feb 20 2022.
%clear all
close all
clc

Is_Save = 1; % 1: generated parameter will be saved
Is_Load = 0; % 1: will abandon generated parameter and load saved file
Is_debugTest = 0; % 1: debug test for first iteration (Is_Load must be 0)

Network_Type = 'small_world';
% random_thread: (connected-ER)linear network with extra links with given probability
% balance_tree: complete binary tree
% abilene: Abilene topology, 11 nodes abilene(DECO version)
% small_world : Watts-Strogatz Small World Graph
% (Deleted)renyi_connected: A Erdos-Renyi graph enhanced to strong connectivity
% GEANT: GEANT network (version by Stratis)
% fog: Fog topology of 4 layers (see DECO)
% LHC: Large Hadron Collider network (see DECO)

eps = 1e-5; % threshold for justify 0
%Network_Type = 'example_1'; % example in figure 3
%% Graph Generating
N_node = 50; % number of nodes

if strcmp(Network_Type,'random_thread')
    p_extralink = 0.07; % graph is a linear network with extra random link
    Adj = Graph_Generator_RandomThread(N_node,p_extralink);
elseif strcmp(Network_Type,'example_1')
    Adj = [0 1 0 1; 1 0 1 1; 0 1 0 1; 1 1 1 0];
elseif strcmp(Network_Type,'balance_tree')
    Adj = Graph_Generator_BalanceTree(N_node);
elseif strcmp(Network_Type,'fog')
    K_server = 2; % number of braches in each layer, totally 4 layers
    K_leaf = 3; % number of leaf at each router
    Adj = Graph_Generator_Fog(K_server,K_leaf);
elseif strcmp(Network_Type,'abilene')
    Adj = Graph_Generator_Abilene();
elseif strcmp(Network_Type,'small_world')
    d_mean = sqrt(N_node); % average degree
    beta = 0.02; % prob of extra link
    Adj = Graph_Generator_SmallWorld(N_node,d_mean,beta);
elseif strcmp(Network_Type,'renyi_connected')
    p_extralink = 0.1; % ER network enhanced to connected
    Adj = Graph_Generator_RenyiConnected(N_node,p_extralink);
elseif strcmp(Network_Type,'GEANT')
    Adj = Graph_Generator_GEANT();
elseif strcmp(Network_Type,'LHC')
    Adj = Graph_Generator_LHC();
end
N_node = length(Adj);
%% applications and tasks
Task_Type = 'random'; % random tasks
N_app = 10;  % number of service chain applications
Ta_min = 2; % mini number of tasks in an application
Ta_max = 2; % mini number of tasks in an application
Source_min = 3;     % minimum active source number
Source_max = 3;   % the maximum input node number of each task
rate_min = 5; % Note: rate_min must greater than eps to avoid unintended task drop
rate_max = 5; % max input rate of task (uniform)
PackSize_sample = [5, 3, 1, 1e-2]; % only works when all Ta = 3
PackSize_sample = PackSize_sample(1:Ta_max+1);

% Ta : row (N_app x 1),
% N_stage: number
Ta = Ta_min -1 + randi(Ta_max +1 -Ta_min, N_app,1);
Da = randperm(N_node,N_app);    % destination
N_stage = sum(Ta) + N_app;
CompWeight = ones(N_node,N_stage);

% Sources: 0-1 maxtrix (N_node x N_app), indecating active sources
% InputRate: matrix (N_node x N_app), each row is input rate for all task of one node
% PackSize: vector (N_stage x 1), order a,k
if strcmp(Task_Type,'random')
    Sources = zeros(N_node, N_app);
    for app = 1:N_app
        source_num = Source_min -1 + randi(Source_max +1 -Source_min);
        source_list = randperm(N_node,source_num);
        for s = source_list
            Sources(s,app) = 1;
        end
    end
    InputRate_init = rate_min + (rate_max - rate_min) * rand(N_node,N_app);
    InputRate = InputRate_init .* Sources;
    PackSize = kron(ones(1,N_app),reshape(PackSize_sample,1,[]));
end

%% phi variables and convert functions
% calculate flow Fij, Gim and Ti+- from routing variable phi, and reverse
% for phi -> f, see function 'Update_Flow' in appendix
% for f -> phi, see function 'FtoPhi' in appendix

%% Delay and Derivative
Delay_type = 'linear'; 
% linear: linear link cost
% queue: queueing delay F/(C-F)

CompCost_type = 'sum_linear'; % 
% sum_linear: sum of all m and linear
% sum_queue: sum of all m and apply queueing delay


Cap_Max_link = 30;
Cap_Min_link = 30;
Delay_para = (Cap_Min_link + (Cap_Max_link-Cap_Min_link) * rand(N_node)) .* Adj; % parameter for delay, e.g. link capacity. Matrix (N x N)
% make the link para symetirc
for node_i = 1:N_node
    for node_j = 1:node_i
        Delay_para(node_i,node_j) = Delay_para(node_j,node_i);
    end
end

% Generate computation capacity according to expoential distribution
Cap_Max_comp = 40;
Cap_Min_comp = 10;
Cap_Exp_para = 15; % the mean value
if strcmp(CompCost_type,'sum_linear')
    CompCost_para = Cap_Min_comp + (Cap_Max_comp-Cap_Min_comp) * rand(N_node,1); % parameter for computation cost, e.g. CPU speed. Column vecror (N x 1)
elseif strcmp(CompCost_type,'sum_queue')
    CompCost_para = Cap_Min_comp + exprnd(Cap_Exp_para,N_node,1);
end
CompCost_para = min(CompCost_para,Cap_Max_comp);

% see function 'Update_Cost' in appendix

%% Save and load topology, parameters
SaveTopoFileName = 'ServiceChainNetwork_SaveTopoPara';
if Is_Save
    save(SaveTopoFileName,'Adj','CompCost_para','CompCost_type','Delay_para','Delay_type','eps',...
        'Source_max','Source_min','Sources','InputRate','Is_debugTest','Is_Load','Is_Save',...
        'N_node','N_app','N_stage','Ta','Da','CompWeight','Network_Type','SaveTopoFileName','PackSize');
end
if Is_Load
    load(SaveTopoFileName);
end
figure(1)
G = graph(Adj);
edges = table2array(G.Edges);
edgeWeight = zeros(size(edges,1),1);
for e = 1:size(edges,1)
    edgeWeight(e) = (Delay_para(edges(e,1),edges(e,2)) + Delay_para(edges(e,2),edges(e,1))) /2 ;%- Cap_Min_link/2; % average link cap for both way
end
p = plot(G,'LineWidth',edgeWeight/5);
p.NodeColor = 'r';
%labelnode(p,1:N_node,'')
for node_i = 1:N_node
    highlight(p,node_i,'MarkerSize',CompCost_para(node_i)/3);
end
saveas(gcf,'NetworkTopology')
saveas(gcf,'NetworkTopology.pdf')

%% Generate initial state
Initial_Type = 'MILP'; % MILP: find a feasibe point using MILP with some arbitrary objective; MILP_RA: random allocation based on a MILP init
Is_Initial_LCOR = 1; % 1: force the initial to be local (or nearest-datasource) computation.

if strcmp(Delay_type, 'queue') % hard upper bound for link flow and compute flow
    LinkCap = Delay_para;
elseif strcmp(Delay_type, 'linear')
    LinkCap = sum(sum(InputRate)) * Delay_para;
end
if strcmp(CompCost_type, 'sum_queue')
    CompCap = CompCost_para;
elseif strcmp(CompCost_type, 'sum_linear')
    CompCap = sum(sum(InputRate)) * CompCost_para;
end

if strcmp(Initial_Type,'MILP') || strcmp(Initial_Type,'MILP_RA')
    Capacity_Saturate_Factor = 0.9; % the maximum allowed fraction of capacity

    %(Adj,N_app,Ta,Da,InputRate,PackSize,Delay_type,Delay_para,CompCost_type,CompCost_para, LinkCap,CompCap,Is_Initial_LCOR)
    [Is_Success,f_init,g_init] = Init_Generator2_MILP(Adj,N_app,Ta,Da,InputRate,PackSize,...
        Delay_type,Delay_para,CompCost_type,CompCost_para,LinkCap*Capacity_Saturate_Factor,CompCap*Capacity_Saturate_Factor,Is_Initial_LCOR);
    if ~Is_Success
        error('ERROR: No Feasible Initial State!  Will abort.\n');
    end
    [Phi_init] = FtoPhi(Adj,N_app,Ta,Da,InputRate,f_init,g_init,eps);
end

Is_Valid = Check_Phi_Valid(Adj,N_app,Ta,Da,InputRate,Phi_init,PackSize,CompWeight,LinkCap, CompCap, eps);
if ~Is_Valid
    error("ERROR: Initial state not valid!");
end
%% Algorithm Parameters
T_MAX = 100; % max iteration number;
EPS_STOP_ratio = 1e-3; % stop if total cost does not decrease more than this ratio

Is_Use_GP = 1;      % if using method: Gradient Projection (scaled matrix = I)
Is_Use_LCOR = 1;   % if using method: LCOR Computation(nearest offloading) + Optimal Routing
Is_Use_LPR = 0;     % if using method: linear programming (relaxed) and rouding, in paper [3]'A Distributed Framework for Task Offloading in Edge Computing Networks of Arbitrary Topology'
Is_Use_SPOO = 1;  % if using method: Shortest path routing (not conjestion aware, use marginal cost at F=0) + optimal offloading


StepSize_GP = 2e-2* ones(1,T_MAX); % use universal stepsize alpha

if Is_Use_GP
    %StepSize_GP = 3e-3 .* (1:T_MAX).^0.6; % use universal stepsize alpha = c/sqrt(t);
    %StepSize_GP = 0.01 .* 1./sqrt(1:T_MAX); % use universal stepsize alpha = c/sqrt(t);
        
    TotalCost_GP = NaN * ones(1,T_MAX);
    Phi_GP_current = Phi_init;
end

if Is_Use_LCOR
    TotalCost_LCOR = NaN * ones(1,T_MAX);
    if Is_Initial_LCOR == 1 % only perform local computation method if start with a local scheme
        Phi_LCOR_current = Phi_init;
    else
        error('ERROR: Not start with LCOR computation!\n')
    end
end
if Is_Use_LPR
    TotalCost_LPR = NaN * ones(1,T_MAX);
end
if Is_Use_SPOO
    TotalCost_SPOO = NaN * ones(1,T_MAX);
    if Is_Initial_LCOR == 1 % note: fixed path should be initialized with local computing and shortest path (marginal at F=0) for result flow
        Phi_SPOO_current = Phi_init;
    else
        error('ERROR: Not start with LCOR computation and shortest path!\n')
    end
end

%% run
Max_Cost = 0; % record of the range of magnitude for plot
Min_Cost = Inf;

time_start = tic;
for iter_t = 1:T_MAX
    if mod(iter_t, 1 ) == 0
        %time_t = toc;
        fprintf('Iteration no. %d, run time %.2f second, ',iter_t,toc(time_start));
    end
    
    % first exam if there is scenario change. If so, re-run some of the methods
    Is_Scenario_Change = 0;
    
    if Is_Use_GP % if using gradient projection algorithm
        % update and save current statues
        [Fij_GP_current,Gi_GP_current,t_GP_current,Is_Loopfree] ...
            = Update_Flow(Adj,N_app,Ta,Da,InputRate,PackSize,CompWeight,Phi_GP_current,eps);
        
        %Fij_GP_current
        %Gi_GP_current
        %t_GP_current

        [Delay_GP_current,Delay_deri_GP_current,CompCost_GP_current,CompCost_deri_GP_current] ...
            = Update_Cost(Delay_type,CompCost_type,Delay_para,CompCost_para,Fij_GP_current,Gi_GP_current,eps);

        TotalCost_GP(iter_t) = sum(sum(Delay_GP_current)) + sum(CompCost_GP_current);

        if mod(iter_t, 1 ) == 0
            %time_t = toc;
            fprintf('total cost = %.3f\n ',TotalCost_GP(iter_t) );
        end

        if TotalCost_GP(iter_t) > Max_Cost
            Max_Cost = TotalCost_GP(iter_t);
        end
        if TotalCost_GP(iter_t) < Min_Cost
            Min_Cost = TotalCost_GP(iter_t);
        end
        [Phi_GP_next] ...
            = Iteration_GP(Adj,N_app,Ta,Da,InputRate,PackSize,CompWeight,StepSize_GP(iter_t),Phi_GP_current,...
            Fij_GP_current,Gi_GP_current,t_GP_current,Delay_deri_GP_current,...
            CompCost_deri_GP_current,Capacity_Saturate_Factor,eps);
        %Adj,N_app,Ta,Da,InputRate,PackSize,CompWeight,StepSize_GP_current,Phi_GP_current,Fij_GP_current,...
    %Gi_GP_current,t_GP_current,Delay_deri_GP_current,CompCost_deri_GP_current,Capacity_Saturate_Factor,eps)
        
        % overwrite
        % Check_Phi_Valid(Adj,N_app,Ta,Da,InputRate,Phi_init,PackSize,CompWeight,LinkCap, CompCap, eps);
        if Check_Phi_Valid(Adj,N_app,Ta,Da,InputRate,Phi_GP_next, PackSize,CompWeight, LinkCap, CompCap, eps)
            Phi_GP_current = Phi_GP_next;
        else
            error('ERROR: unvalid phi in method: GP.\n');
        end        
               
    end
    
    
    if Is_Use_LCOR % if using local (neareast to data source) computation
        % update and save current statues
        [Fij_LCOR_current,Gi_LCOR_current,t_LCOR_current,Is_Loopfree] ...
            = Update_Flow(Adj,N_app,Ta,Da,InputRate,PackSize,CompWeight,Phi_LCOR_current,eps);
        %Adj,N_app,Ta,Da,InputRate,PackSize,CompWeight,Phi,eps
        
        [Delay_LCOR_current,Delay_deri_LCOR_current,CompCost_LCOR_current,CompCost_deri_LCOR_current] ...
            = Update_Cost(Delay_type,CompCost_type,Delay_para,CompCost_para,Fij_LCOR_current,Gi_LCOR_current,eps);
        %Delay_type,CompCost_type,Delay_para,CompCost_para,Fij,Gi,eps

        TotalCost_LCOR(iter_t) = sum(sum(Delay_LCOR_current)) + sum(CompCost_LCOR_current);
        
        if TotalCost_LCOR(iter_t) > Max_Cost
            Max_Cost = TotalCost_LCOR(iter_t);
        end
        if TotalCost_LCOR(iter_t) < Min_Cost
            Min_Cost = TotalCost_LCOR(iter_t);
        end
        
        % update variable: using SGP for phi+, but keep phi-
        [Phi_LCOR_next] ...
            = Iteration_LCOR(Adj,N_app,Ta,Da,InputRate,PackSize,CompWeight,StepSize_GP(iter_t),Phi_LCOR_current,...
            Fij_LCOR_current,Gi_LCOR_current,t_LCOR_current,Delay_deri_LCOR_current,...
            CompCost_deri_LCOR_current,Capacity_Saturate_Factor,eps);
        %(Adj,N_app,Ta,Da,InputRate,PackSize,CompWeight,StepSize_LCOR_current,Phi_LCOR_current,Fij_LCOR_current,Gi_LCOR_current,...
            %t_LCOR_current,Delay_deri_LCOR_current,CompCost_deri_LCOR_current,Capacity_Saturate_Factor,eps)
        
        % overwrite
        %Check_Phi_Valid(Adj,N_app,Ta,Da,InputRate,Phi,PackSize,CompWeight,LinkCap, CompCap, eps)
        if Check_Phi_Valid(Adj,N_app,Ta,Da,InputRate,Phi_LCOR_next,PackSize,CompWeight,LinkCap,CompCap, eps)
            Phi_LCOR_current = Phi_LCOR_next;
        else
            error('ERROR: unvalid phi in method: LCOR.\n');
        end
        
        
    end
    
    if Is_Use_LPR % if using linear programming and rounding provided in eq(11) in paper [3]
        % first compare current topology and task with previous, only perform LPR if t=1 or there's a topology change
        
        if Is_Scenario_Change == 1 % if any of Adj,Task and InputRate has changed, reperform LPR
            [Phi_minus_LPR_current,Phi_plus_LPR_current] = Update_Phi_LPR(Adj,M,Task,InputRate,a_m,b_overhead,...
                Delay_type,Delay_para,CompCost_type,CompCost_para,LinkCap,CompCap,Capacity_Saturate_Factor,eps);
            
            if ~Check_Phi_Valid(Adj,M,Task,InputRate,a_m,b_overhead,Phi_minus_LPR_current,Phi_plus_LPR_current, LinkCap, CompCap, eps)
                error('ERROR: unvalid phi in method: LPR.\n');
            end
            
            [LinkFlow_data_LPR_current,LinkFlow_result_LPR_current,CompFlow_LPR_current,t_minus_LPR_current,t_plus_LPR_current,Is_Loopfree] ...
                = Update_Flow(Adj,M,Task,a_m,b_overhead,InputRate,Phi_minus_LPR_current,Phi_plus_LPR_current,eps);
            LinkFlow_LPR_current = LinkFlow_data_LPR_current + LinkFlow_result_LPR_current;
            [Delay_LPR_current,Delay_deri_LPR_current,CompCost_LPR_current,CompCost_deri_LPR_current] ...
                = Update_Cost(Delay_type,CompCost_type,Delay_para,CompCost_para,LinkFlow_LPR_current,CompFlow_LPR_current,eps);
            
            TotalCost_LPR(iter_t) = sum(sum(Delay_LPR_current)) + sum(sum(CompCost_LPR_current));
            
            Adj_previous = Adj;
            Task_previous = Task;
            InputRate_previous = InputRate;
            
        else % if non of these changed, keep previous result
            TotalCost_LPR(iter_t) = TotalCost_LPR(iter_t-1);
        end
        
        if TotalCost_LPR(iter_t) > Max_Cost
            Max_Cost = TotalCost_LPR(iter_t);
        end
        if TotalCost_LPR(iter_t) < Min_Cost
            Min_Cost = TotalCost_LPR(iter_t);
        end
    end
    
    if Is_Use_SPOO % if using shortest path + optimal offloading
        % Note: the initial state is already set to shortest path, so the updating would use SGP, but keep the blocked nodes
        % to the nodes where initially no flow
        % Note: the new initial state and blocked nodes will be computed if the scenario has changes
        
              
        % first calculate current total cost
        [Fij_SPOO_current,Gi_SPOO_current,t_SPOO_current,Is_Loopfree] ...
            = Update_Flow(Adj,N_app,Ta,Da,InputRate,PackSize,CompWeight,Phi_SPOO_current,eps);
        
         [Delay_SPOO_current,Delay_deri_SPOO_current,CompCost_SPOO_current,CompCost_deri_SPOO_current] ...
            = Update_Cost(Delay_type,CompCost_type,Delay_para,CompCost_para,Fij_SPOO_current,Gi_SPOO_current,eps);

        TotalCost_SPOO(iter_t) = sum(sum(Delay_SPOO_current)) + sum(CompCost_SPOO_current);
        if TotalCost_SPOO(iter_t) > Max_Cost
            Max_Cost = TotalCost_SPOO(iter_t);
        end
        if TotalCost_SPOO(iter_t) < Min_Cost
            Min_Cost = TotalCost_SPOO(iter_t);
        end
        
        Is_Blocked_SPOO = Block_Keep_Path(Adj,N_app,Ta,Phi_SPOO_current,eps);

        % then update phi according to given block matrix
        [Phi_SPOO_next] ...
            = Iteration_SPOO(Adj,N_app,Ta,Da,InputRate,PackSize,CompWeight,StepSize_GP(iter_t),Phi_SPOO_current,...
            Fij_SPOO_current,Gi_SPOO_current,t_SPOO_current,Delay_deri_SPOO_current,...
            CompCost_deri_SPOO_current,Capacity_Saturate_Factor,Is_Blocked_SPOO,eps);
        
        % update phi
        if Check_Phi_Valid(Adj,N_app,Ta,Da,InputRate,Phi_SPOO_next,PackSize,CompWeight,LinkCap,CompCap, eps)
            Phi_SPOO_current = Phi_SPOO_next;
        else
            error('ERROR: unvalid phi in method: SPOO.\n');
        end
    end
    
end

%% plot
figure(2)
if Is_Use_GP
    plot(TotalCost_GP(1:end),'b-','DisplayName','Gradient Prpjection','LineWidth',1.5)
    TotalCost_GP_Final = TotalCost_GP(end)
    %[Lifetime_data_GP,Lifetime_result_GP] = Average_Liftime(Fij_GP_current,Delay_GP_current,Task,InputRate,a_m);
    %Liftime_ratio_GP = Lifetime_data_GP/Lifetime_result_GP;
    %Liftime_and_ratio_GP = [Lifetime_data_GP,Lifetime_result_GP,Liftime_ratio_GP];
    hold on
end
if Is_Use_LPR
    plot(TotalCost_LPR(1:end),'m-','DisplayName','Linear Program Rounded [3]','LineWidth',1.5)
    TotalCost_LPR_Final = TotalCost_LPR(end)
    [Lifetime_data_LPR,Lifetime_result_LPR] = Average_Liftime(LinkFlow_data_LPR_current,LinkFlow_result_LPR_current,Delay_LPR_current,Task,InputRate,a_m);
    Liftime_ratio_LPR = Lifetime_data_LPR/Lifetime_result_LPR;
    Liftime_and_ratio_LPR = [Lifetime_data_LPR,Lifetime_result_LPR,Liftime_ratio_LPR];
    hold on
end
if Is_Use_LCOR
    plot(TotalCost_LCOR(1:end),'g-','DisplayName','LCOR Computation Optimal Routing','LineWidth',1.5)
    TotalCost_LCOR_Final = TotalCost_LCOR(end)
    %[Lifetime_data_LCOR,Lifetime_result_LCOR] = Average_Liftime(LinkFlow_data_LCOR_current,LinkFlow_result_LCOR_current,Delay_LCOR_current,Task,InputRate,a_m);
    %Liftime_ratio_LCOR = Lifetime_data_LCOR/Lifetime_result_LCOR;
    %Liftime_and_ratio_LCOR = [Lifetime_data_LCOR,Lifetime_result_LCOR,Liftime_ratio_LCOR];
    hold on
end
if Is_Use_SPOO
    plot(TotalCost_SPOO(1:end),'k-','DisplayName','Shortest Path Optimal Offloading','LineWidth',1.5)
    TotalCost_SPOO_Final = TotalCost_SPOO(end)
    %[Lifetime_data_SPOO,Lifetime_result_SPOO] = Average_Liftime(LinkFlow_data_SPOO_current,LinkFlow_result_SPOO_current,Delay_SPOO_current,Task,InputRate,a_m);
    %Liftime_ratio_SPOO = Lifetime_data_SPOO/Lifetime_result_SPOO;
    %Liftime_and_ratio_SPOO = [Lifetime_data_SPOO,Lifetime_result_SPOO,Liftime_ratio_SPOO];
    hold on
end

xlabel('Iteration');
ylabel('Total Cost');
axis([1 T_MAX Min_Cost-(Max_Cost-Min_Cost)*0.2 Max_Cost+(Max_Cost-Min_Cost)*0.2])
legend
hold off
saveas(gcf,'Output')
saveas(gcf,'Output.pdf')

%% Save data and analyze
SaveDataFileName = 'RouteComputeNetwork_Data';
save(SaveDataFileName);


%% Appendix: Functions
function [Phi_GP_next] = ...
    Iteration_GP(Adj,N_app,Ta,Da,InputRate,PackSize,CompWeight,StepSize_GP_current,Phi_GP_current,Fij_GP_current,...
    Gi_GP_current,t_GP_current,Delay_deri_GP_current,CompCost_deri_GP_current,Capacity_Saturate_Factor,eps)
% Main iteration for Gradient Projection variable updating

% update marginal delay pD/pr and pD/pt via broadcasting
[Margin_node_GP_current,Margin_link_GP_current] ...
    = Broadcast_Margin(Adj,N_app,Ta,Da,PackSize,CompWeight,Phi_GP_current,...
    Delay_deri_GP_current,CompCost_deri_GP_current,eps);
%Margin_node_GP_current
%Adj,N_app,Ta,Da,PackSize,CompWeight,Phi,Delay_deri,CompCost_deri,eps

% calculate bloacked nodes
[Is_Blocked_GP_current] ...
    = Update_Blocked(Adj,N_app,Ta,Phi_GP_current,Margin_node_GP_current,eps);
%Is_Blocked_GP_current
%Adj,N_app,Ta,Phi,Margin_node,eps

% new phi after projection
[Phi_GP_next] ...
    = Update_Phi_GP(Adj,N_app,Ta,Da,Phi_GP_current,Margin_link_GP_current,Is_Blocked_GP_current,StepSize_GP_current,eps);
%Update_Phi_GP(Adj,N_app,Ta,Da,Phi,Margin_link,Is_Blocked,stepsize,eps)
end

function [Phi_LCOR_next] = ...
    Iteration_LCOR(Adj,N_app,Ta,Da,InputRate,PackSize,CompWeight,stepsize,Phi_LCOR_current,Fij_LCOR_current,Gi_LCOR_current,...
    t_LCOR_current,Delay_deri_LCOR_current,CompCost_deri_LCOR_current,Capacity_Saturate_Factor,eps)
% Main iteration for local (nearest to data source) computation
% use SGP method to the result flow, while keep the data flow unchanged, since the initiator is made to local computation

% update marginal delay pD/pr and pD/pt via broadcasting
[Margin_node_LCOR_current,Margin_link_LCOR_current] ...
    = Broadcast_Margin(Adj,N_app,Ta,Da,PackSize,CompWeight,Phi_LCOR_current,...
    Delay_deri_LCOR_current,CompCost_deri_LCOR_current,eps);
%Adj,N_app,Ta,Da,PackSize,CompWeight,Phi,Delay_deri,CompCost_deri,eps

% calculate bloacked nodes
[Is_Blocked_LCOR_current] ...
    = Update_Blocked(Adj,N_app,Ta,Phi_LCOR_current,Margin_node_LCOR_current,eps);
%Adj,N_app,Ta,Phi,Margin_node,eps

% new phi after projection
[Phi_LCOR_next] ...
    = Update_Phi_LCOR(Adj,N_app,Ta,Da,Phi_LCOR_current,Margin_link_LCOR_current, ...
    Is_Blocked_LCOR_current,stepsize,eps);
%Adj,N_app,Ta,Da,Phi,Margin_link,Is_Blocked,stepsize,eps
end

function [Phi_SPOO_next] = ...
    Iteration_SPOO(Adj,N_app,Ta,Da,InputRate,PackSize,CompWeight,stepsize,Phi_SPOO_current,Fij_SPOO_current,Gi_SPOO_current,...
    t_SPOO_current,Delay_deri_SPOO_current,CompCost_deri_SPOO_current,Capacity_Saturate_Factor,Is_Blocked_SPOO,eps)
% shortest path optimal offloading
[Margin_node_SPOO_current,Margin_link_SPOO_current] ...
    = Broadcast_Margin(Adj,N_app,Ta,Da,PackSize,CompWeight,Phi_SPOO_current,...
    Delay_deri_SPOO_current,CompCost_deri_SPOO_current,eps);
%Adj,N_app,Ta,Da,PackSize,CompWeight,Phi,Delay_deri,CompCost_deri,eps

% new phi after projection
% note : SPOO is GP confined to shortest path, so only change the blocked set
[Phi_SPOO_next] ...
    = Update_Phi_GP(Adj,N_app,Ta,Da,Phi_SPOO_current,Margin_link_SPOO_current, ...
    Is_Blocked_SPOO,stepsize,eps);
%Adj,N_app,Ta,Da,Phi,Margin_link,Is_Blocked,stepsize,eps
end

function [Phi_minus,Phi_plus] = ...
    Update_Phi_LPR(Adj,M,Task,InputRate,a_m,b_overhead, Delay_type,Delay_para,CompCost_type,CompCost_para, LinkCap,CompCap,Capacity_Saturate_Factor,eps)
% implementation of centrailized version of paper [3]

% Step 1: use linear programming to solve the relaxed L-JoSRAT, (11)(12) in paper [3]
% note that the algorithm do not consider result flow, so DataFlow_Saturate_Factor should be loose (e.g. 0.8)
MAX_PATH_NUM = 3;
DataFlow_Saturate_Factor = Capacity_Saturate_Factor - 0.2; % upper bound of ratio of capacity for data flow on links (only applies on link flow, not computation flow)
[x_bar_LPR,Subtask,N_subtask,paths_LPR] = LJOSRAT_LPR(Adj,Task,InputRate,Delay_type,Delay_para,CompCost_type,CompCost_para,...
    LinkCap,CompCap,Capacity_Saturate_Factor,DataFlow_Saturate_Factor,MAX_PATH_NUM,eps);

% Step 2: round the fractional solution to integer using Alg2 in [3]
x_LPR = Rounding_LPR(x_bar_LPR,paths_LPR,Adj,Subtask,N_subtask,LinkCap,CompCap,Capacity_Saturate_Factor,MAX_PATH_NUM,eps);

% Step 3: transform the integer solution to x variable, assgin the result flow using min-hop, and project to feasible set
[Phi_minus,Phi_plus] = xToFlow_LPR(x_LPR,paths_LPR,Adj,M,Task,InputRate,Subtask,N_subtask,...
    a_m,b_overhead,LinkCap,CompCap,Capacity_Saturate_Factor,MAX_PATH_NUM,eps);
end

function [x_bar_LPR, Subtask, N_subtask, paths_LPR] = ...
    LJOSRAT_LPR(Adj,Task,InputRate, Delay_type,Delay_para,CompCost_type,CompCost_para, LinkCap,CompCap,Capacity_Saturate_Factor,DataFlow_Saturate_Factor,MAX_PATH_NUM,eps)
% solve the linear programming problem (11)(12) in [3].
% note: do not consider RAM resource; each node could be server;
% tasks are split to single-source-single-dest subtasks, and each subtask could only be offload to one server;
% for each subtask, pick the first MAX_PATH_NUM shortest path (hop number) as candidate
% use Capacity * DataFlow_Saturate_Factor for upper bound of link flow, and Capacity * Capacity_Saturate_Factor for computation flow
% x_bar_LPR: column vector of length N_subtask * N_node * MAX_PATH_NUM
% Subtask: matrix representing subtasks split from task: (N_subtask x 5), each row consist of (source, destination, computation type, input rate,task id)
% paths_LPR: matrix of (N_subtask*N_node*MAX_PATH_NUM x N_node), stroing all paths
fprintf('LJOSRAT_LPR: begin solving LPR...\n');
N_node = length(Adj);
N_task = size(Task,1);
%N_subtask = sum(sum(InputRate > eps));
% split task into subtasks
Subtask = zeros(N_task*N_node,5); % list of subtasks, each row consist of (source, destination, computation type, input rate, task id)
subt_index = 1;
for t_index = 1:N_task
    t_dest = Task(t_index,1);
    t_m = Task(t_index,2);
    for node_n = 1:N_node
        if InputRate(node_n,t_index) > eps % only count if input rate is larger that eps
            Subtask(subt_index,1) = node_n;
            Subtask(subt_index,2) = t_dest;
            Subtask(subt_index,3) = t_m;
            Subtask(subt_index,4) = InputRate(node_n,t_index);
            Subtask(subt_index,5) = t_index;
            subt_index = subt_index +1;
        end
    end
end
N_subtask = subt_index -1;
Subtask = sortrows(Subtask,4,'descend'); % sort the tasks in the descend order of input rate


% set LP parameters: min f*x, s.t. Ax <= b
Length_x_bar = N_subtask * N_node * MAX_PATH_NUM; % number of x varibables, equal to number of paths
lb_LPR = zeros(1,Length_x_bar);
ub_LPR = ones(1,Length_x_bar);

N_cons_LinkCap = N_node * N_node;
N_cons_CompCap = N_node;
N_cons_Offload = N_subtask;
N_cons_LinkBig = N_node * N_node;
N_cons_CompBig = N_node;
N_cons = N_cons_LinkCap + N_cons_CompCap + N_cons_Offload + N_cons_LinkBig + N_cons_CompBig;
A_LPR = zeros(N_cons , Length_x_bar);
b_LPR = zeros(N_cons,1);

% Generate shortest paths: using package given by 'https://www.mathworks.com/matlabcentral/fileexchange/32513-k-shortest-path-yen-s-algorithm'
N_paths = N_subtask * N_node * MAX_PATH_NUM; % ordered as #of subtask, # of node as computation server(containing source itself but not activated), and # of path from source to server
paths_LPR = zeros(N_paths,N_node); % a row is a path, end with 0; if a row start with 0, it is not an active path.
Adj_modified = zeros(size(Adj)); % adj matrix with inf at no-link cases
for node_i = 1:N_node
    for node_j = 1:N_node
        if Adj(node_i,node_j) > eps
            Adj_modified(node_i,node_j) = Adj(node_i,node_j);
        else
            Adj_modified(node_i,node_j) = Inf;
        end
    end
end
paths_LPR_t_struct = zeros(N_node * MAX_PATH_NUM ,N_node,N_subtask);
parfor t_index = 1:N_subtask
    fprintf('Preparing paths for subtask %d\n',t_index);
    t_source = Subtask(t_index,1);
    paths_LPR_t = zeros( N_node * MAX_PATH_NUM ,N_node);
    for node_server = 1:N_node
        [shortestPaths, totalCosts] = kShortestPath(Adj_modified, t_source, node_server, MAX_PATH_NUM); % find the shortest paths w.r.t. hop numbers
        for path_id = 1:length(shortestPaths)
            %path_pos = (t_index-1)*N_node*MAX_PATH_NUM + (node_server-1)*MAX_PATH_NUM + path_id;
            path_pos_t = (node_server-1)*MAX_PATH_NUM + path_id;
            path_temp = cell2mat(shortestPaths(path_id));
            %paths_LPR(path_pos,1:length(path_temp)) = path_temp;
            paths_LPR_t(path_pos_t,1:length(path_temp)) = path_temp;
        end
    end
    paths_LPR_t_struct(:,:,t_index) = paths_LPR_t;
end
for t_index = 1:N_subtask
    paths_LPR((t_index-1)*N_node*MAX_PATH_NUM + 1: t_index*N_node*MAX_PATH_NUM , : ) = paths_LPR_t_struct(:,:,t_index);
end

% first the objective f: f_tip = cost(local compute t) - cost(offload t to i via path p)
% note: the link delay is computed by only transmitting task t, no other traffic
f_LPR = zeros(1,Length_x_bar);
for t_index = 1:N_subtask
    % compute the w_\inf (comp cost for local computation)
    %rate_comp = Subtask(t_index,4);
    rate_comp = Subtask(t_index,4) * 1.5; % lift the computation load to consider congestion
    %rate_link = Subtask(t_index,4);
    rate_link = Subtask(t_index,4) / 1.5; % supress link flow to encourage offloading
    para = CompCost_para(Subtask(t_index,1));
    if strcmp(CompCost_type,'sum_queue')
        if rate_comp/para < Capacity_Saturate_Factor % if could computed locally
            CompCost_LCOR = rate_comp / (para - rate_comp);
        else % if could NOT computed locally, assign a very high cost
            CompCost_LCOR = 1/eps;
        end
    elseif strcmp(CompCost_type,'sum_linear')
        CompCost_LCOR = rate_comp * para;
    else
        error('ERROR: undefined computation type!\n');
    end
    % compute and assign entries for all paths
    for node_server = [1:(Subtask(t_index,1)-1) (Subtask(t_index,1)+1):N_node] % all nodes but not the requester
        % compute cost at server
        para = CompCost_para(node_server);
        if strcmp(CompCost_type,'sum_queue')
            if rate_comp/para < Capacity_Saturate_Factor % if could computed locally
                CompCost_Server = rate_comp / (para - rate_comp);
            else % if could NOT computed locally, assign a very high cost
                CompCost_Server = 1/eps;
            end
        elseif strcmp(CompCost_type,'sum_linear')
            CompCost_Server = rate_comp * para;
        else
            error('ERROR: undefined computation type!\n');
        end
        
        for path_index = 1:MAX_PATH_NUM % calculate link costs for data flow from source to server
            path_pos = (t_index-1)*N_node*MAX_PATH_NUM + (node_server-1)*MAX_PATH_NUM + path_index;
            path = paths_LPR(path_pos,:);
            path_len = length(find(path));
            if path_len >= 2 % only consider active paths
                DelaySum = 0;
                for link_index = 1:path_len-1
                    link_i = path(link_index);
                    link_j = path(link_index +1);
                    para = Delay_para(link_i,link_j);
                    
                    if strcmp(Delay_type,'queue')
                        if rate_link/para < Capacity_Saturate_Factor % if could computed locally
                            DelaySum = DelaySum + rate_link / (para - rate_link);
                        else % if could NOT computed locally, assign a very high cost
                            DelaySum = DelaySum + 1/eps;
                        end
                    elseif strcmp(Delay_type,'linear')
                        DelaySum = DelaySum + rate_link * para;
                    else
                        error('ERROR: undefined computation type!\n');
                    end
                end
                Offload_Gain = CompCost_LCOR - (CompCost_Server + DelaySum);
                x_bar_pos = path_pos;
                f_LPR(x_bar_pos) = -1 * Offload_Gain; % maximizing gain is minimizing -gain
            end
        end
    end
end

% then the constraint matrix A and b
% note: 4 kinds of constraints: link capacity, computation capacity, offloading at most 1, and the 'big task' constraints; ordered as follows

% link cap
% \sum_{p:ij \in p} r_p x_p <= cij;
% for all p, find its r_p and assgin to correct entry
for t_index = 1:N_subtask
    rate_comp = Subtask(t_index,4);
    t_source = Subtask(t_index,1);
    for node_server = [1:(t_source-1) (t_source+1):N_node] % for all possible server
        for path_index = 1:MAX_PATH_NUM % for all possible paths
            path_pos = (t_index-1)*N_node*MAX_PATH_NUM + (node_server-1)*MAX_PATH_NUM + path_index;
            path = paths_LPR(path_pos,:);
            path_len = length(find(path));
            if path_len >= 2 % only consider active paths
                x_bar_pos = path_pos;
                % link caps
                for link_index = 1:path_len-1 % for all links on path
                    link_i = path(link_index);
                    link_j = path(link_index +1);
                    Cons_LinkCap_pos = (link_i-1)*N_node + link_j; % # of linkcap ij in constraints
                    LinkCap_ij = LinkCap(link_i,link_j);
                    A_LPR(Cons_LinkCap_pos,x_bar_pos) = rate_comp;
                    b_LPR(Cons_LinkCap_pos) = LinkCap_ij * DataFlow_Saturate_Factor; % would be assigned multiple times, but not matter
                end
            end
        end
    end
end
% comp cap
% \sum_{p:p_end = i} r_p x_p <= ci;
for t_index = 1:N_subtask
    rate_comp = Subtask(t_index,4);
    t_source = Subtask(t_index,1);
    for node_server = [1:(t_source-1) (t_source+1):N_node]
        Cons_CompCap_pos = N_cons_LinkCap + node_server;
        b_LPR(Cons_CompCap_pos) = CompCap(node_server) * Capacity_Saturate_Factor;
        for path_index = 1:MAX_PATH_NUM % for all possible paths
            path_pos = (t_index-1)*N_node*MAX_PATH_NUM + (node_server-1)*MAX_PATH_NUM + path_index;
            path = paths_LPR(path_pos,:);
            path_len = length(find(path));
            if path_len >= 2 % only consider active paths
                x_bar_pos = path_pos;
                A_LPR(Cons_CompCap_pos,x_bar_pos) = rate_comp;
            end
        end
    end
end
% offloading sum up to 1
% \sum_{p:p \in t} x_p <= 1;
for t_index = 1:N_subtask
    Cons_Offload_pos = N_cons_LinkCap + N_cons_CompCap + t_index;
    b_LPR(Cons_Offload_pos) = 1;
    t_source = Subtask(t_index,1);
    for node_server = [1:(t_source-1) (t_source+1):N_node]
        for path_index = 1:MAX_PATH_NUM % for all possible paths
            path_pos = (t_index-1)*N_node*MAX_PATH_NUM + (node_server-1)*MAX_PATH_NUM + path_index;
            path = paths_LPR(path_pos,:);
            path_len = length(find(path));
            if path_len >= 2 % only consider active paths
                x_bar_pos = path_pos;
                A_LPR(Cons_Offload_pos,x_bar_pos) = 1;
            end
        end
    end
end
% big task for links
% \sum_{p: r_p > cij/2} x_p <= 1, forall ij
for t_index = 1:N_subtask
    rate_comp = Subtask(t_index,4);
    t_source = Subtask(t_index,1);
    for node_server = [1:(t_source-1) (t_source+1):N_node] % for all possible server
        for path_index = 1:MAX_PATH_NUM % for all possible paths
            path_pos = (t_index-1)*N_node*MAX_PATH_NUM + (node_server-1)*MAX_PATH_NUM + path_index;
            path = paths_LPR(path_pos,:);
            path_len = length(find(path));
            if path_len >= 2 % only consider active paths
                x_bar_pos = path_pos;
                % link caps
                for link_index = 1:path_len-1 % for all links on path
                    link_i = path(link_index);
                    link_j = path(link_index +1);
                    Cons_LinkBig_pos = N_cons_LinkCap + N_cons_CompCap + N_cons_Offload + (link_i-1)*N_node + link_j; % # of linkcap ij in constraints
                    LinkCap_ij = LinkCap(link_i,link_j);
                    if rate_comp >= LinkCap_ij * DataFlow_Saturate_Factor /2 % if t is a big task wrt ij
                        A_LPR(Cons_LinkBig_pos,x_bar_pos) = 1;
                        b_LPR(Cons_LinkBig_pos) = 1;
                    end
                end
            end
        end
    end
end
% big task for servers
% \sum_{p: r_p > ci/2} x_p <=1, forall i
for t_index = 1:N_subtask
    rate_comp = Subtask(t_index,4);
    t_source = Subtask(t_index,1);
    for node_server = [1:(t_source-1) (t_source+1):N_node] % for all possible server
        Cons_CompBig_pos = N_cons_LinkCap + N_cons_CompCap + N_cons_Offload + N_cons_LinkBig + node_server;
        b_LPR(Cons_CompBig_pos) = 1;
        CompCap_server = CompCap(node_server);
        if rate_comp >= CompCap_server * Capacity_Saturate_Factor / 2 % if the current t takes more that 1/2 at server, it is a big task
            for path_index = 1:MAX_PATH_NUM % for all possible paths
                path_pos = (t_index-1)*N_node*MAX_PATH_NUM + (node_server-1)*MAX_PATH_NUM + path_index;
                path = paths_LPR(path_pos,:);
                path_len = length(find(path));
                if path_len >= 2 % only consider active paths
                    x_bar_pos = path_pos;
                    A_LPR(Cons_CompBig_pos,x_bar_pos) = 1;
                end
            end
        end
    end
end

% finally, carry out the LP
[x_bar_LPR,fval,exitflag,output] = linprog(f_LPR,A_LPR,b_LPR,[],[],lb_LPR,ub_LPR);
fprintf('LJOSRAT_LPR: LPR finished\n');
end

function x_int_LPR = ...
    Rounding_LPR(x_bar_LPR,paths_LPR,Adj,Subtask,N_subtask,LinkCap,CompCap,Capacity_Saturate_Factor,MAX_PATH_NUM,eps)
% round the fractional solution x_bar_LPR into integer, using algorithm 2 in paper [3].
N_node = length(Adj);
%N_subtask = size(Subtask,1);
% Step 1: for each task, pick the largest p to offload. If no x_p is non-zero, then compute at data source
x_int_LPR = zeros(size(x_bar_LPR));
for t_index = 1:N_subtask
    x_sum = 0; % accumulate the sum of x_p for subtask t
    x_max = 0;
    x_max_pos = 0; % the path p with max x_p for subtask t
    t_source = Subtask(t_index,1);
    for node_server = [1:(t_source-1) (t_source+1):N_node]
        for path_index = 1:MAX_PATH_NUM
            path_pos = (t_index-1)*N_node*MAX_PATH_NUM + (node_server-1)*MAX_PATH_NUM + path_index;
            if path_pos > length(paths_LPR)
                error('!')
            end
            path = paths_LPR(path_pos,:);
            path_len = length(find(path));
            if path_len >= 2 % only consider active paths
                x_pos = path_pos;
                x_sum = x_sum + x_bar_LPR(x_pos);
                if x_bar_LPR(x_pos) > x_max
                    x_max = x_bar_LPR(x_pos);
                    x_max_pos = x_pos;
                end
            end
        end
    end
    if x_sum < 0.5 % if the sum of x_p is small, decide to compute locally
        continue
    else % if sum of x_p is large, offload to the maximum x_p
        x_int_LPR(x_max_pos) = 1;
    end
end

% Step 2: Alternation phase and final phase
% Omit here

end

function [Phi_minus,Phi_plus] = ...
    xToFlow_LPR(x_LPR,paths_LPR,Adj,M,Task,InputRate,Subtask,N_subtask,a_m,b,LinkCap,CompCap,Capacity_Saturate_Factor,MAX_PATH_NUM,eps)
% calculate the link flows according to given routing-offloading variable x_LPR
% Assign the data flow accroding to x_LPR, and the result flow the min-hop path
% After the assginment above, project f-,f+ and g to the feasible set (using QP)
% After projection, converte the flow variable into phi variable, for convenience of further SGP iteration on the result flow (if needed)
N_node = length(Adj);
N_task = size(Task,1);
fprintf('xToFlow_LPR: start\n');
% Step 1: construct the f- and g flow
% note: the order of flow f is i-j-t, for g is i-t, and for x_LPR is t-v-p
length_f_minus = N_node * N_node * N_task;
f_minus_flow = zeros(length_f_minus,1);

length_g = N_node *  N_task;
g_flow = zeros(length_g,1);

parfor subt_index = 1:N_subtask
    t_index = Subtask(subt_index,5);
    subt_source = Subtask(subt_index,1);
    rate = Subtask(subt_index,4);
    subt_server = Subtask(subt_index,1); % actual server for t, defualt is local computation
    f_minus_flow_t = zeros(length_f_minus,1);
    g_flow_t =  zeros(length_g,1);
    for node_server = [1:(subt_source -1) (subt_source+1):N_node]
        for path_index = 1:MAX_PATH_NUM
            path_pos = (subt_index-1)*N_node*MAX_PATH_NUM + (node_server-1)*MAX_PATH_NUM + path_index;
            x_LPR_pos = path_pos;
            if abs( x_LPR(x_LPR_pos) -1) <= eps % for all p such that x_p = 1
                path = paths_LPR(path_pos,:);
                path_len = length(find(path));
                for link_index = 1:path_len-1 % for all links on path p, add rate to f-ijt
                    link_i = path(link_index);
                    link_j = path(link_index +1);
                    f_minus_pos = (link_i-1)*N_node*N_task + (link_j-1)*N_task + t_index;
                    %f_minus_flow(f_minus_pos) = f_minus_flow(f_minus_pos) + rate;
                    f_minus_flow_t(f_minus_pos) = f_minus_flow_t(f_minus_pos) + rate;
                end
                % save the actual server
                subt_server = node_server;
            end
        end
    end
    % for the corr. server, add rate to git
    g_pos = (subt_server -1)*N_task + t_index;
    %g_flow(g_pos) = g_flow(g_pos) + rate;
    g_flow_t(g_pos) = g_flow_t(g_pos) + rate;
    
    f_minus_flow = f_minus_flow + f_minus_flow_t;
    g_flow = g_flow + g_flow_t;
end

% Step 2: calculate the min-hop feasible result flow using LP on f_plus_flow (since computation is fixed, LP still valid even if b neq 0)
% min fx , s.t. Aleq*x = bleq, Aeq*x = beq

% f_plus_flow = zeros(length_f_plus,1);
% length_x_flow = length_f_minus + length_f_plus + length_g;
% x_flow = zeros(1,length_x_flow);
% length_f_plus = length_f_minus;

G = digraph(Adj);
 [sOut,tOut] = findedge( G );
Edge = [sOut tOut]; % list of all edges (directed), each row is the [node_i,node_j]
N_edge = size(Edge,1);
length_f_plus = N_edge * N_task;

f_LP = ones(1,length_f_plus); % to obtain min-hop, i.e. minimizes sum of flow
length_cons_eq = N_node * N_task; % flow conservation of f+ for each i,t
Aeq_LP = zeros(length_cons_eq,length_f_plus);
beq_LP = zeros(length_cons_eq,1);
for t_index = 1:N_task
    t_dest = Task(t_index,1);
    t_m = Task(t_index,2);
    for node_i = 1:N_node
        % \sum_j fij  - \sum_j fji = (a*g+b)  (if i is not dest)
        % \sum_j fij = 0 (if i is dest)
        cons_eq_pos = (node_i-1)*N_task + t_index;
        
        % first calculate inject flow due to local computation
        g_pos = (node_i -1)*N_task + t_index;
        InjectFlow = g_flow(g_pos) * a_m(t_m);
        if g_flow(g_pos) >= eps
            InjectFlow = InjectFlow + b; % inject flow = a*g + b
        end
        
        % then for all in-going and out-going flow
        for node_j = find(Adj(node_i,:)) % for j: ij is link
            %fij_pos = (node_i-1)*N_node*N_task + (node_j-1)*N_task + t_index;
            %fij_pos = (node_i-1)*N_node*N_task + (node_j-1)*N_task + t_index;
            
            [tf, edge_index]=ismember([node_i node_j],Edge,'rows');
            fij_pos = (edge_index-1)*N_task + t_index;
            
            Aeq_LP(cons_eq_pos,fij_pos) = 1;
        end
        if node_i == t_dest % if is the dest
            beq_LP(cons_eq_pos) = 0;
        else % if is not dest
            beq_LP(cons_eq_pos) = InjectFlow;
            for node_j = reshape(find(Adj(:,node_i)),1,[]) % for j: ji is link
                %fji_pos = (node_j-1)*N_node*N_task + (node_i-1)*N_task + t_index;
                
                 [tf, edge_index]=ismember([node_j node_i],Edge,'rows');
                fji_pos = (edge_index-1)*N_task + t_index;
                
                Aeq_LP(cons_eq_pos,fji_pos) = -1;
            end
        end
    end
end

%length_cons_leq = N_node * N_node; % link capacity consrtaints for each link
length_cons_leq = N_edge; % link capacity consrtaints for each link

Aleq_LP = zeros(length_cons_leq,length_f_plus);
bleq_LP = zeros(length_cons_leq,1);
for node_i = 1:N_node
    for node_j = find(Adj(node_i,:)) % for all ij is link
        % \sum_t f+ijt <= cij - \sum_t f-ij
        
        %cons_leq_pos = (node_i-1)*N_node + node_j;
        
         [tf, edge_index]=ismember([node_i node_j],Edge,'rows');
        cons_leq_pos = edge_index;
        
        bleq_LP(cons_leq_pos) = LinkCap(node_i,node_j) * Capacity_Saturate_Factor; % start with capacity, will be subtracted with f-ij
        for t_index = 1:N_task
            %fijt_plus_pos = (node_i-1)*N_node*N_task + (node_j-1)*N_task + t_index;
            %fijt_minus_pos = fijt_plus_pos
            fijt_plus_pos = (edge_index -1)*N_task + t_index;       
            fijt_minus_pos = (node_i-1)*N_node*N_task + (node_j-1)*N_task + t_index;
            
            Aleq_LP(cons_leq_pos,fijt_plus_pos) = 1;
            bleq_LP(cons_leq_pos) = bleq_LP(cons_leq_pos) - f_minus_flow(fijt_minus_pos);
        end
    end
end
lb_LP = zeros(1,length_f_plus);
ub_LP = [];
% carry out LP to find min-hop feasible f+
[f_plus_edge,fval,exitflag,output] = linprog(f_LP,Aleq_LP,bleq_LP,Aeq_LP,beq_LP,lb_LP,ub_LP);

length_f_plus_output = N_node^2*N_task;
 f_plus_flow = zeros(length_f_plus_output,1);
    %g_comp_init = zeros(length_g_output,1);
    for t_index = 1:N_task
        for edge_index = 1:N_edge
            node_i = Edge(edge_index,1);
            node_j = Edge(edge_index,2);
            pos_f_edge = (edge_index -1)*N_task + t_index;
            pos_f_flow = (node_i -1)*N_node*N_task + (node_j-1)*N_task + t_index;
            f_plus_flow(pos_f_flow) = f_plus_edge(pos_f_edge);
        end
    end

% Step 3: convert to phi
[Phi_minus,Phi_plus] =  FtoPhi(Adj,Task,InputRate,a_m,b,f_minus_flow,f_plus_flow,g_flow,eps);

Is_Valid = Check_Phi_Valid(Adj,M,Task,InputRate,a_m,b,Phi_minus,Phi_plus, LinkCap, CompCap, eps);
if  ~Is_Valid
    error('ERROR: Unfeasible solution by LPR!\n');
end
fprintf('LPR: Done converting to phi\n');
end

function [Is_Blocked] = ...
    Block_Keep_Path(Adj,N_app,Ta,Phi,eps)
% Calculate the set of blocked nodes for given phi, to keep the routing path
% (namely, block all links that do not sit on shorest (min-hop) path to the destination)
% Note: local computation is never blocked.
% Is_Blocked_minus, Is_Blocked_plus: column vector (N*N*N_task x 1), (i,j,t) denotes if j is blocked to i wrt task t.
N_node = length(Adj);
N_stage = sum(Ta) + N_app;
Sa = Ta + 1; % number of stages for each app

length_block_vec = N_node * N_node *N_stage;
Is_Blocked = ones(length_block_vec,1);

for node_i = 1:N_node
    for node_j = find(Adj(node_i,:))
        % for all links ij, note that the un-link ij is already blocked by initializing to 1.
        for app = 1:N_app
            for k = 0:Ta(app)
                s_index = sum(Sa(1:app-1))+k +1;
                phi_pos = (node_i-1)*(N_node+1)*N_stage + node_j*N_stage + s_index;
                block_vec_pos = (node_i-1)*N_node*N_stage + (node_j-1)*N_stage + s_index;
                if Phi(phi_pos) >= eps % if there is non-zeros phi, unblock
                    Is_Blocked(block_vec_pos) = 0;
                end
            end
        end
    end
end

end


function [Scale_Mat_minus,Scale_Mat_plus]  = ...
    Update_Scale_Matrix(Adj,Task,StepSize_GP,Phi_minus,Phi_plus,Is_Blocked_minus,Is_Blocked_plus, t_minus,t_plus,...
    Margin_link_minus,Margin_link_plus,Scale_Type,Capacity_Saturate_Factor,eps)
% Calculate scaling matrix M- and M+ for all i,(d,m)
% Note: Since M- and M+ are diagnal matrix, we only store diagnal elements
% Scale_Mat_minus: column vector (N*(1+N)*N_task x 1), ordered as i,j,t, where {j} is the diag elements of M_i(t).
% Scale_Mat_plus: column vector (N*N*N_task x 1), ordered as i,j,t, where {j} is the diag elements of M_i(t).
% StepSize_GP_current: scalar, universal stepsize of current time slot
N_node = length(Adj);
N_task = size(Task,1);

length_Scale_Mat_minus = N_node*(1+N_node)*N_task;
length_Scale_Mat_plus = N_node*N_node*N_task;
Scale_Mat_minus = zeros(length_Scale_Mat_minus,1);
Scale_Mat_plus = zeros(length_Scale_Mat_plus,1);
if strcmp(Scale_Type,'GP')
    % if using simple gradient descent, M = tidm/ alphaidm * diag{1,..0,..1}, only 0 is on pos j that minimized delta_ijdm.
    % see eq. (37) in Yufang paper
    for t_index = 1:N_task
        for node_i = 1:N_node
            % calculate diagonal elements for M-idm and M+idm
            
            % first minus
            % find the minimum marginal out-link
            min_margin_minus_pos = -1;
            min_margin_minus = Inf;
            for node_j = 0:N_node % note: start from 0
                Margin_minus_pos = (node_i -1)*(1+N_node)*N_task + node_j*N_task + t_index; % pos of deltaij
                if Margin_link_minus(Margin_minus_pos) < min_margin_minus
                    min_margin_minus = Margin_link_minus(Margin_minus_pos);
                    min_margin_minus_pos = Margin_minus_pos;
                end
            end
            if (min_margin_minus_pos < 0) || (min_margin_minus == Inf)
                error('ERROR: Fail to find minimum link marginal');
            end
            % assgin entries
            for node_j = 0:N_node
                M_minus_pos = (node_i -1)*(1+N_node)*N_task + node_j*N_task + t_index;
                if M_minus_pos == min_margin_minus_pos
                    Scale_Mat_minus(M_minus_pos) = 0;
                else
                    t_minus_pos = (node_i-1)*N_task + t_index;
                    Scale_Mat_minus(M_minus_pos) = t_minus(t_minus_pos) / StepSize_GP;
                end
            end
            
            % then plus
            min_margin_plus_pos = -1;
            min_margin_plus = Inf;
            for node_j = 1:N_node % note: start from 1
                Margin_plus_pos = (node_i -1)*N_node*N_task + (node_j-1)*N_task + t_index; % pos of deltaij
                if Margin_link_plus(Margin_plus_pos) < min_margin_plus
                    min_margin_plus = Margin_link_plus(Margin_plus_pos);
                    min_margin_plus_pos = Margin_plus_pos;
                end
            end
            if (min_margin_plus_pos < 0) || (min_margin_plus == Inf)
                error('ERROR: Fail to find minimum link marginal');
            end
            % assgin entries
            for node_j = 1:N_node
                M_plus_pos = (node_i -1)*N_node*N_task + (node_j-1)*N_task + t_index;
                if M_plus_pos == min_margin_plus_pos
                    Scale_Mat_plus(M_plus_pos) = 0;
                else
                    t_plus_pos = (node_i-1)*N_task + t_index;
                    Scale_Mat_plus(M_plus_pos) = t_plus(t_plus_pos) / StepSize_GP;
                end
            end
            
        end
    end
elseif strcmp(Scale_Type,'SGP')
    % for scaled gradient projection, M-, M+ is given in eq.(41) Yufang paper, or see page 11 in overleaf
    % M+ = t+/2 * diag{Aij(D0) + |O(i)\B+idm|*h+jdm*A(D0)}_j in O(i)\B+idm
    % Aij(D0) = max_{Dij < D0} D''ij.
    % For queueing delay D = F/(C-F), D'' = 2C/(C-F).
    % If at initial state we assume Fij <= aij*Cij, then D'' <= 2/(1-a).  aij: Capacity_Saturate_Factor
    % A(D0) = max_{ij} Aij(D0)
    % hjdm-/+: max hop path length from g to destination/sink
    
    % calculate Aij(D0) and A(D0)
    Aij_D0_plus = zeros(N_node,N_node);
    Aij_D0_minus = zeros(N_node,N_node+1);
    A_D0 = 0;
    Link_Saturate_Factor = min(Capacity_Saturate_Factor,0.9);
    for node_i = 1:N_node
        for node_j = find(Adj(node_i,:))
            Aij_D0_plus(node_i,node_j) = 2/(1-Link_Saturate_Factor);
            if Aij_D0_plus(node_i,node_j) > A_D0
                A_D0 = Aij_D0_plus(node_i,node_j);
            end
        end
    end
    for node_i = 1:N_node
        for node_j = [0 find(Adj(node_i,:))]
            Aij_D0_minus(node_i,node_j+1) = 2/(1-Link_Saturate_Factor);
            if Aij_D0_minus(node_i,node_j+1) > A_D0
                A_D0 = Aij_D0_plus(node_i,node_j+1);
            end
        end
    end
    
    for t_index = 1:N_task
        
        % calculate hidm-/+ using bfs
        % reconstruct adj matrix for positive phi
        Adj_plus_t = zeros(N_node,N_node);
        Adj_minus_t = zeros(N_node,N_node);
        for node_i = 1:N_node
            for node_j = find(Adj(node_i,:))
                Phi_plus_pos = (node_i-1)*N_node*N_task + (node_j-1)*N_task + t_index;
                Phi_minus_pos = (node_i-1)*(1+N_node)*N_task + node_j*N_task + t_index;
                Adj_plus_t(node_i,node_j) = ( Phi_plus(Phi_plus_pos) > eps);
                Adj_minus_t(node_i,node_j) = ( Phi_minus(Phi_minus_pos) > eps);
            end
        end
        h_minus = ones(1,N_node); % for minus flow, start from 1 to include local computation
        h_plus = zeros(1,N_node);
        % first plus
        Is_traversed = zeros(1,N_node);
        Is_traverse_next = zeros(1,N_node); % mark nodes that not visited but all out-neighbors visited
        t_dest = Task(t_index,1);
        Is_traverse_next(t_dest) = 1;
        while sum(Is_traversed) < N_node
            for node_i = find(Is_traverse_next) % add max hop by 1
                for node_j = find(Adj_plus_t(node_i,:))
                    h_plus(node_i) = max(h_plus(node_i),h_plus(node_j)+1);
                end
                Is_traversed(node_i) = 1;
            end
            Is_traverse_next = zeros(1,N_node);
            for node_i = find(Is_traversed == 0) % find to-visit node
                Is_traverse_next(node_i) = (min(Is_traversed(find(Adj_plus_t(node_i,:)))) == 1 );
                % if all i's out-neighbor has been traversed, but i is not
            end
        end
        %then minus, start from nodes only forward to local cpu
        Is_traversed = zeros(1,N_node);
        Is_traverse_next = zeros(1,N_node);
        t_sink = find(sum(Adj_minus_t,2) == 0);
        Is_traverse_next(t_sink) = 1;
        while sum(Is_traversed) < N_node
            for node_i = find(Is_traverse_next) % add max hop by 1
                for node_j = find(Adj_minus_t(node_i,:))
                    h_minus(node_i) = max(h_minus(node_i),h_minus(node_j)+1);
                end
                Is_traversed(node_i) = 1;
            end
            Is_traverse_next = zeros(1,N_node);
            for node_i = find(Is_traversed == 0) % find to-visit node
                Is_traverse_next(node_i) = (min(Is_traversed(find(Adj_minus_t(node_i,:)))) == 1 );
                % if all i's out-neighbor has been traversed, but i is not
            end
        end
        %h_minus
        %h_plus
        
        % update matrix entries
        for node_i = 1:N_node
            % first update M+
            t_plus_pos = (node_i -1)*N_task + t_index;
            Is_blocked_plus_pos_list = ((node_i -1)*N_node*N_task + t_index) :N_task ...
                : ((node_i -1)*N_node*N_task + (N_node-1)*N_task + t_index); % list of Bijdm+ for all j
            Unblocked_size_plus = N_node - sum(Is_Blocked_plus(Is_blocked_plus_pos_list)); % size of O(i)\Bidm+
            for node_j = find(Adj(node_i,:))
                Is_blocked_plus_pos = (node_i -1)*N_node*N_task + (node_j-1)*N_task + t_index;
                M_plus_pos = (node_i -1)*N_node*N_task + (node_j-1)*N_task + t_index;
                if Is_Blocked_plus(Is_blocked_plus_pos) == 0 % if j is in O(i)\Bidm+
                    %Scale_Mat_plus(M_plus_pos) = t_plus(t_plus_pos) / 2 * ( Aij_D0_plus(node_i,node_j) +  A_D0 * Unblocked_size_plus * h_plus(node_j));
                    Scale_Mat_plus(M_plus_pos) = t_plus(t_plus_pos) / 2 * ( 0*Aij_D0_plus(node_i,node_j) +  2 * Unblocked_size_plus * h_plus(node_j));
                end
            end
            % then update M-, note: include j = 0
            t_minus_pos = (node_i -1)*N_task + t_index;
            Is_blocked_minus_pos_list = ((node_i -1)*N_node*N_task + t_index) :N_task ...
                : ((node_i -1)*N_node*N_task + (N_node-1)*N_task +t_index); % note: i0 is not in blocked matrix
            Unblocked_size_minus = 1+N_node - sum(Is_Blocked_minus(Is_blocked_minus_pos_list)); % size of O(i)\Bidm-
            for node_j = [0 find(Adj(node_i,:))]
                
                M_minus_pos = (node_i-1)*(1+N_node)*N_task + node_j*N_task + t_index;
                if node_j == 0
                    Scale_Mat_minus(M_minus_pos) = t_minus(t_minus_pos) / 2 ...
                        * ( Aij_D0_minus(node_i,node_j+1) +  A_D0 * Unblocked_size_minus * 0);
                else
                    Is_blocked_minus_pos = (node_i -1)*N_node*N_task + (node_j-1)*N_task + t_index;
                    if Is_Blocked_minus(Is_blocked_minus_pos) == 0 % if j is in O(i)\Bidm-
                        %Scale_Mat_minus(M_minus_pos) = t_minus(t_minus_pos) / 2 * ( Aij_D0_minus(node_i,node_j+1) +  A_D0 * Unblocked_size_minus * h_minus(node_j));
                        Scale_Mat_minus(M_minus_pos) = t_minus(t_minus_pos) / 2 * ( 0*Aij_D0_minus(node_i,node_j+1) +  2 * Unblocked_size_minus * h_minus(node_j));
                    end
                end
                
            end
        end
    end
else
    error('ERROR: Wrong Scale Type');
end

end

function [Margin_node,Margin_link] = Broadcast_Margin(Adj,N_app,Ta,Da,PackSize,CompWeight,Phi,Delay_deri,CompCost_deri,eps)
% Calculate the marginal costs pD/pt, simulating broadcast
% Phi: column vecto (N*(N+1)*N_stage x 1)
% Delay_deri, matrix (N x N)
% Compcostderi, Vector (N x 1)
% Margin_node: matrix (N x N_stage)
% Margin_link: column vector ( N*(N+1)*N_task x 1 )
% coeff_t = 0.1;
% offset_t = max(sum(InputRate,2));

N_node = length(Adj);
N_stage = sum(Ta) + N_app;
Sa = Ta + 1; % number of stages for each app

length_Phi = N_node * (N_node + 1) * N_stage;

% broadcast start from final stage
Margin_node = zeros(N_node,N_stage);

for app = 1:N_app
    da = Da(app);
    for k = Ta(app):-1:0
        % for each stage, using BFS-like starting from destination
        s_index = sum(Sa(1:app-1))+k +1;

        Is_traversed = zeros(1,N_node); % Marks for completed nodes
        Is_taverse_next = zeros(1,N_node); % Marks for nodes that will determin pD/pt next iteration
        DAG_stage = zeros(N_node,N_node); % Adj matrix for directed graph of current stage

        % construncting DAG
        for node_i = 1:N_node
            for node_j = 1:N_node
                Phi_pos = (node_i -1)* (N_node+1) * N_stage + node_j * N_stage + s_index;
                if Phi(Phi_pos) > eps % mark the link with phi > 0
                    if Adj(node_i,node_j) == 0
                        error('ERROR: Positive phi+ on non-link (%d,%d)!\n',node_i,node_j)
                    end
                    DAG_stage(node_i,node_j) = 1;
                end
            end
        end

        % finding starting point
        % if k = Ta, start with da
        % if k neq Ta, start with nodes with no out-neighbor
        Is_taverse_next = zeros(1,N_node);
        if k == Ta(app)
            Is_taverse_next(da) = 1;
        else
            for node_i = 1:N_node
                if sum(DAG_stage(node_i,:)) == 0
                    Is_taverse_next(node_i) = 1;
                end
            end
        end

        % traverse nodes with all out-neighbor been traversed
        while sum(Is_traversed) < N_node % continue if there is node left
            if sum(Is_taverse_next) == 0 % if not finished but no nodes to traverse next
                error('ERROR: Broadcast abortion!\n')
            end
            %Is_taverse_nextrf
            for node_i = find(Is_taverse_next) % for all nodes to traverse
                for node_j = find(DAG_stage(node_i,:))  % all out neighbor of i
                    % update pD/pt+i and delta_ij+

                    if Is_traversed(node_j) == 0
                        error('ERROR: Traversing node is not ready!\n');
                    end
                    Phi_pos = (node_i -1)* (N_node+1) * N_stage + node_j * N_stage + s_index;

                    % pD/pt_i = \sum_j phiij+ (L*D'ij + pD/pt+j)
                    Margin_node(node_i,s_index) = Margin_node(node_i,s_index) + ...
                        Phi(Phi_pos) * (Margin_node(node_j,s_index) + PackSize(s_index) * Delay_deri(node_i,node_j));

                end

                % if not the final stage, also count for comput marginal w*Ci'
                if k ~= Ta(app)
                    margin_next_stage = Margin_node(node_i,s_index + 1);
                    Phi_i0_pos = (node_i -1)* (N_node+1) * N_stage + s_index;
                    Margin_node(node_i,s_index) = Margin_node(node_i,s_index) + ...
                        Phi(Phi_i0_pos) * (margin_next_stage + CompWeight(node_i,s_index) * CompCost_deri(node_i));
                end

                Is_traversed(node_i) = 1;
            end

            % then mark all nodes to be visit next iter: nodes not traversed but with out-neighbors all traversed.
            Is_taverse_next = zeros(1,N_node);
            for node_i = find(Is_traversed == 0)
                Is_will_traverse = 1;
                for node_j = find(DAG_stage(node_i,:)) % for all j: (i,j) in DAG
                    if Is_traversed(node_j) == 0
                        Is_will_traverse = 0;
                    end
                end
                if Is_will_traverse == 1
                    Is_taverse_next(node_i) = 1;
                end
            end
        end
    end
end

% finally globally compute link/compute marginals
Margin_link = zeros(length_Phi,1);
for app = 1:N_app
    for k = 0:Ta(app)
        s_index = sum(Sa(1:app-1))+k +1;
        for node_i = 1:N_node
            if k~=Ta(app)
                delta_i0_pos = (node_i-1)*(1+N_node)*N_stage + s_index;
                Margin_link(delta_i0_pos) = ... % delta_i0 = w*C'im +  pD/pt(k+1)
                    CompWeight(node_i,s_index) * CompCost_deri(node_i) + Margin_node(node_i,s_index+1);
            end
            for node_j = find(Adj(node_i,:)) % for all j: (i,j) is link, delta_ij = L*D'ij + pDpT(j)
                delta_pos = (node_i-1)*(1+N_node)*N_stage + node_j*N_stage + s_index;
                Margin_link(delta_pos) = PackSize(s_index) * Delay_deri(node_i,node_j) + Margin_node(node_j,s_index);% delta_ij- = D'ij + pD/prj
            end
        end
    end
end
end

function [Phi_new] = Update_Phi_GP(Adj,N_app,Ta,Da,Phi,Margin_link,Is_Blocked,stepsize,eps)
% 1-order method, no scaling matrix
% Phi_minus_new: column vector (N*(N+1)*N_task x 1);
% Phi_plus_new: column vector (N*N*N_task x 1);

N_node = length(Adj);
N_stage = sum(Ta) + N_app;
Sa = Ta + 1; % number of stages for each app
Phi_new = zeros(size(Phi));

for app = 1:N_app
    da = Da(app);
    for k = 0:Ta(app)
        s_index = sum(Sa(1:app-1))+k +1;
        for node_i = 1:N_node
            % eq (8) (9) (10)
            % first, find min delta_ij and identify N_i
            % then, calculate e_ij, 
            % then, calculate Delta_phi_ij for j not min, sum up to S
            % finally, devide S to min-delta-directions

             %fprintf('Examing node_i = %d\n',node_i);
            if k == Ta(app) && node_i == da
                continue
            end

            % extract phi, delta, block
            Phi_iak = zeros(1,N_node+1);
            Phi_new_iak = Phi_iak;
            Margin_link_iak = zeros(1,N_node+1); % \delta_ij(a,k) for j in [0 V]
            Is_blocked_iak = zeros(1,N_node+1); % note: first item is comp
            for node_j = 0:N_node
                marg_pos = (node_i -1)*(N_node +1)*N_stage + node_j*N_stage + s_index;
                Margin_link_iak(node_j+1) = Margin_link(marg_pos);
                Phi_iak(node_j+1) = Phi(marg_pos);
            end
            sum_phi_old = sum(Phi_iak);
            if k == Ta(app) % if final stage, block to comp
                Is_blocked_iak(1) = 1;
            else
                Is_blocked_iak(1) = 0;
            end
            for node_j = 1:N_node
                block_pos = (node_i -1)*N_node*N_stage + (node_j -1)*N_stage + s_index;
                Is_blocked_iak(node_j+1) = Is_Blocked(block_pos);
            end

            % find min delta, Ni
            min_delta = min(Margin_link_iak(Is_blocked_iak ~= 1));
            if isempty(min_delta)
                return
            end
            Diff = abs(Margin_link_iak - min_delta);
            N_iak = find( Diff < eps ); % note: start with j=0

            % calculate e_ij, for blocked e = phi,; for block else e = phi-min;
            e_iak = zeros(1,N_node+1);
            for node_j = 0:N_node
                if Is_blocked_iak(node_j +1) == 1
                    e_iak(node_j+1) = 0;
                else
                    e_iak(node_j+1) = Margin_link_iak(node_j+1) - min_delta;
                end
            end

            % update phi for non-min directions
            S_iak = 0;
            for node_j = 0:N_node
                
                if isempty(find(N_iak == node_j+1)) % if j is not in N_i
                    %fprintf('Non-min-node %d\n',node_j);
                    if Is_blocked_iak(node_j+1) == 1 % if blocked, delet all phi
                        Phi_new_iak(node_j+1) = 0;
                          %fprintf("phi_pre = %f, phi_cur = %f, d = %f.\n",Phi_iak(node_j+1),Phi_new_iak(node_j+1),Phi_iak(node_j+1))
                        S_iak = S_iak + Phi_iak(node_j+1);
                    else % if not blocked, substract min{phi, alpha * e}
                        d_ijak = min(Phi_iak(node_j+1), stepsize*e_iak(node_j+1));
                        Phi_new_iak(node_j+1) = Phi_iak(node_j+1) - d_ijak;
                        S_iak = S_iak + d_ijak;
                        %fprintf("phi_pre = %f, phi_cur = %f, d = %f.\n",Phi_iak(node_j+1),Phi_new_iak(node_j+1),d_ijak)
                    end
                end
            end
            % update phi for min directions
            Num_min = length(N_iak);
            inc = S_iak/Num_min;
            for j_pos = N_iak
                %fprintf('Min-node %d\n',j_pos-1);
                Phi_new_iak(j_pos) = Phi_iak(j_pos) + inc;
            end

            if min(Phi_new_iak) < 0
                return
            end
            
            sum_phi_new = sum(Phi_new_iak);
            if abs(sum_phi_new - sum_phi_old) > eps
                error("Sum of phi has changed!\n");
            end
            
            if node_i == 5 && s_index == 1
                %fprintf('min_delta = %.3f\n',min_delta)
                %Margin_link_iak
                %Is_blocked_iak
                %Phi_iak
                %Phi_new_iak
            end

            % assign entries
            for node_j = 0:N_node
                Phi_pos = (node_i -1)*(1+N_node)*N_stage + node_j*N_stage + s_index;
                Phi_new(Phi_pos) = Phi_new_iak(node_j+1);
            end
        end
    end
end

end

function [Phi_new] = Update_Phi_LCOR(Adj,N_app,Ta,Da,Phi,Margin_link,Is_Blocked,stepsize,eps)
% Similar to function 'Update_Phi_GP', but keep phi for stage 0 and update phi+
%options =  optimoptions(@quadprog,'Display','off');


N_node = length(Adj);
N_stage = sum(Ta) + N_app;
Sa = Ta + 1; % number of stages for each app
Phi_new = zeros(size(Phi));

for app = 1:N_app
    da = Da(app);

    % first copy stage 0 phi to new vect
    s_index = sum(Sa(1:app-1))+1;
    for node_i = 1:N_node
        sum_phi = 0;
        phi_ia0 = zeros(1,N_node+1);
        for node_j = 0:N_node
            Phi_pos = (node_i -1)*(1+N_node)*N_stage + node_j*N_stage + s_index;
            phi_ia0(node_j+1) = Phi(Phi_pos);
            Phi_new(Phi_pos) = phi_ia0(node_j+1);
            sum_phi = sum_phi + Phi_new(Phi_pos);
        end
        if sum_phi~= 1
            
        end
    end

    for k = 1:Ta(app)
        s_index = sum(Sa(1:app-1))+k +1;
        for node_i = 1:N_node
            % eq (8) (9) (10)
            % first, find min delta_ij and identify N_i
            % then, calculate e_ij, 
            % then, calculate Delta_phi_ij for j not min, sum up to S
            % finally, devide S to min-delta-directions

             %fprintf('Examing node_i = %d\n',node_i);
            if k == Ta(app) && node_i == da
                continue
            end

            % extract phi, delta, block
            Phi_iak = zeros(1,N_node+1);
            Phi_new_iak = Phi_iak;
            Margin_link_iak = zeros(1,N_node+1); % \delta_ij(a,k) for j in [0 V]
            Is_blocked_iak = zeros(1,N_node+1); % note: first item is comp
            for node_j = 0:N_node
                marg_pos = (node_i -1)*(N_node +1)*N_stage + node_j*N_stage + s_index;
                Margin_link_iak(node_j+1) = Margin_link(marg_pos);
                Phi_iak(node_j+1) = Phi(marg_pos);
            end
            sum_phi_old = sum(Phi_iak);
            if k == Ta(app) % if final stage, block to comp
                Is_blocked_iak(1) = 1;
            else
                Is_blocked_iak(1) = 0;
            end
            for node_j = 1:N_node
                block_pos = (node_i -1)*N_node*N_stage + (node_j -1)*N_stage + s_index;
                Is_blocked_iak(node_j+1) = Is_Blocked(block_pos);
            end

            % find min delta, Ni
            min_delta = min(Margin_link_iak(Is_blocked_iak ~= 1));
            if isempty(min_delta)
                return
            end
            Diff = abs(Margin_link_iak - min_delta);
            N_iak = find( Diff < eps ); % note: start with j=0

            % calculate e_ij, for blocked e = phi,; for block else e = phi-min;
            e_iak = zeros(1,N_node+1);
            for node_j = 0:N_node
                if Is_blocked_iak(node_j +1) == 1
                    e_iak(node_j+1) = 0;
                else
                    e_iak(node_j+1) = Margin_link_iak(node_j+1) - min_delta;
                end
            end

            % update phi for non-min directions
            S_iak = 0;
            for node_j = 0:N_node
                
                if isempty(find(N_iak == node_j+1)) % if j is not in N_i
                    %fprintf('Non-min-node %d\n',node_j);
                    if Is_blocked_iak(node_j+1) == 1 % if blocked, delet all phi
                        Phi_new_iak(node_j+1) = 0;
                          %fprintf("phi_pre = %f, phi_cur = %f, d = %f.\n",Phi_iak(node_j+1),Phi_new_iak(node_j+1),Phi_iak(node_j+1))
                        S_iak = S_iak + Phi_iak(node_j+1);
                    else % if not blocked, substract min{phi, alpha * e}
                        d_ijak = min(Phi_iak(node_j+1), stepsize*e_iak(node_j+1));
                        Phi_new_iak(node_j+1) = Phi_iak(node_j+1) - d_ijak;
                        S_iak = S_iak + d_ijak;
                        %fprintf("phi_pre = %f, phi_cur = %f, d = %f.\n",Phi_iak(node_j+1),Phi_new_iak(node_j+1),d_ijak)
                    end
                end
            end
            % update phi for min directions
            Num_min = length(N_iak);
            inc = S_iak/Num_min;
            for j_pos = N_iak
                %fprintf('Min-node %d\n',j_pos-1);
                Phi_new_iak(j_pos) = Phi_iak(j_pos) + inc;
            end

            if min(Phi_new_iak) < 0
                return
            end
            
            sum_phi_new = sum(Phi_new_iak);
            if abs(sum_phi_new - sum_phi_old) > eps
                error("Sum of phi has changed!\n");
            end

            % assign entries
            for node_j = 0:N_node
                Phi_pos = (node_i -1)*(1+N_node)*N_stage + node_j*N_stage + s_index;
                Phi_new(Phi_pos) = Phi_new_iak(node_j+1);
            end
        end
    end
end
end


function [Is_Blocked] = Update_Blocked(Adj,N_app,Ta,Phi,Margin_node,eps)
% Calculate the blocked nodes accroding to pD/pt.
% Is_Blocked: (N*N*N_stage x 1), (i,j,t) denotes if j is blocked to i wrt task t.
% Method:   First calculate improper nodes by BFS (i.e., containing an improper path to destination/sink)
%           Then assign entries to blocked_mat: either no ij link, or phiij = 0 and j is improper.
N_node = length(Adj);
N_stage = sum(Ta) + N_app;
Sa = Ta + 1; % number of stages for each app

length_Is_Blocked = N_node * N_node * N_stage;
Is_Blocked = zeros(length_Is_Blocked,1);
for app = 1:N_app % separate for each task
    for k = 0:Ta(app)
        s_index = sum(Sa(1:app-1))+k +1;

        Is_improper = zeros(1,N_node);
        Is_improper_link = zeros(N_node,N_node); % improper links indicator
        Phi_s = zeros(N_node,N_node);
        % first mark all improper links
        for node_i = 1:N_node
            for node_j = find(Adj(node_i,:)) % all links i,j
                Phi_pos = (node_i-1)*(1+N_node)*N_stage + node_j*N_stage + s_index;
                Phi_s(node_i,node_j) = Phi(Phi_pos);
                if (Phi_s(node_i,node_j) > eps) ... % first plus: if phiij+ >0 and pD/pt+i < pD/pt+j
                        && (Margin_node(node_j,s_index) - Margin_node(node_i,s_index) > eps)
                    Is_improper_link(node_i,node_j) = 1;
                end               
            end
        end

        % then BSF (flood for at most N iterations)
        for iter = 1:N_node
            Is_changed = 0;
            % plus
            for node_i = find(Is_improper == 0) % exam all nodes not marked as improper,
                for node_j = find(Adj(node_i,:))
                    if Is_improper_link(node_i,node_j) == 1 % mark if it has improper out-link
                        %[node_i,node_j]
                        Is_improper(node_i) = 1;
                        Is_changed = 1;
                        break;
                    elseif Is_improper(node_j) == 1% mark if it has improper out-neighbor
                        Is_improper(node_i) = 1;
                        Is_changed = 1;
                        break;
                    end
                end
            end
            %Is_improper
            if Is_changed == 0
                break;
            end
        end

        % then assign to output
        for node_i = 1:N_node
            for node_j = 1:N_node
                Is_Blocked_pos = (node_i-1)*N_node*N_stage + (node_j-1)*N_stage + s_index;
                if Adj(node_i,node_j) == 0 % if not a link, block both - and +
                    Is_Blocked(Is_Blocked_pos) = 1;
                else % if ij is a link
                    Phi_pos = (node_i -1)*(1+N_node)*N_stage + node_j*N_stage + s_index;
                    if Phi(Phi_pos) < eps % only block nodes with no current flow
                        if Margin_node(node_j,s_index) - Margin_node(node_i,s_index) > eps % if there is link ij, but j has higher marginal
                            Is_Blocked(Is_Blocked_pos) = 1;
                        end
                        if Is_improper(node_j) == 1
                            Is_Blocked(Is_Blocked_pos) = 1;
                        end
                    end
                end
            end
        end
        if  0
            s_index
            Phi_s
            Margin_node
            Is_improper
            Is_improper_link
        end
    end

end
end

function [Task,InputRate] = ...
    Task_Generator_random(Adj,N_app,T_a,rate_min,rate_max,inpute_nodes_number_max,eps)
% generates application destination, packet size, etc.
N_node = length(Adj);
Task = zeros(N_task,2);
for i = 1:N_task
    Task(i,1) = randi(N_node); % each row of Task is a (d,m)
    Task(i,2) = randi(M);
end

% random input rate
input_nodes_number = randi(inpute_nodes_number_max,1,N_task); % input node number of each task
InputRate = zeros(N_node,N_task);  % order: i,t
for i_task = 1:N_task
    input_n_vec = randperm(N_node, input_nodes_number(i_task));
    for i_node = 1:length(input_n_vec)
        InputRate(input_n_vec(i_node),i_task) = rate_min + rand() * (rate_max-rate_min);
    end
end
end

function Adj = ...
    Graph_Generator_RandomThread(N_node,p_extralink)
%if Network_Type == 'random_thread' % random network grown from a linear topology, each link with a prob p_extralink
Adj = (rand(N_node) < p_extralink);
for i = 1:N_node-1
    Adj(i,i+1) = 1;
end
for i = 1:N_node
    Adj(i,i) = 0;
    for j = 1:i-1
        Adj(i,j) = Adj(j,i);
    end
end
%end
end

function Adj = ...
    Graph_Generator_BalanceTree(N_node)
% generate a complete binary tree
% round the input N_node to nearest 2^K
K = round(log2(N_node));

Adj = zeros(2^K-1);
if K == 1
    return;
end
% link: every node i is connected to 2i and 2i+1, if i is not leaf
for node_i = 1:2^(K-1)-1
    Adj(node_i,2*node_i) = 1;
    Adj(node_i,2*node_i+1) = 1;
end
Adj = Adj + Adj.'; % make symmetric
end

function Adj = ...
    Graph_Generator_Fog(K_server,K_leaf)
% generate a complete 4-depth K-tree, and linear connection on the same layer (without leaf)
N_node = 1 + K_server + K_server^2 + K_leaf*(K_server^2);
Adj = zeros(N_node);
% first layer  branches
for node_j = 1+(1:K_server)
    Adj(1,node_j) = 1;
end
% first layer cross link
for node_i = 2:K_server
    Adj(node_i,node_i+1) = 1;
end
%Adj(2,1+K_server) = 1;
% Second layer branches
for branch_i = 1:K_server
    node_i = 1 + branch_i;
   for branch_j = 1:K_server
       node_j = 1 + K_server + (branch_i-1)*K_server + branch_j;
       Adj(node_i,node_j) = 1;
   end
end
% Second layer cross links
for node_i = 2+K_server:K_server + K_server^2
    Adj(node_i,node_i+1) = 1;
end
%Adj(2+K_server,1+K_server+K_server^2) = 1;
% Third layer (leaf)
for branch_i = 1:K_server^2
    node_i = (1+K_server) + branch_i;
    for branch_j = 1:K_leaf
        node_j = 1 + K_server + K_server^2 + (branch_i-1)*K_leaf + branch_j;
        Adj(node_i,node_j) = 1;
    end
end
% leaf cross links
for branch_i = 1:K_server^2
    for leaf = 1+K_server+K_server^2+(branch_i-1)*K_leaf+1:1+K_server+K_server^2+(branch_i-1)*K_leaf+K_leaf-1
        Adj(leaf,leaf+1) = 1;
    end
end


Adj = Adj + Adj.'; % make symmetric
end

function Adj = ...
    Graph_Generator_Abilene()
% generate the Abilene Network (DECO version)
Adj = zeros(11);

Adj(1,2) = 1;
Adj(1,3) = 1;
Adj(2,3) = 1;
Adj(2,4) = 1;
Adj(4,5) = 1;
Adj(5,6) = 1;
Adj(3,6) = 1;
Adj(5,7) = 1;
Adj(6,8) = 1;
Adj(7,8) = 1;
Adj(8,9) = 1;
Adj(7,10) = 1;
Adj(10,11) = 1;
Adj(9,11) = 1;

Adj = Adj + Adj.'; % make symmetric
end

function Adj = ...
    Graph_Generator_SmallWorld(N_node,d_mean,beta)
% Watts-Strogatz Small world , see https://www.mathworks.com/help/matlab/math/build-watts-strogatz-small-world-graph-model.html
N = N_node;
K = floor(d_mean /2);
s = repelem((1:N)',1,K);
t = s + repmat(1:K,N,1);
t = mod(t-1,N)+1;

% Rewire the target node of each edge with probability beta
for source=1:N    
    switchEdge = rand(K, 1) < beta;
    
    newTargets = rand(N, 1);
    newTargets(source) = 0;
    newTargets(s(t==source)) = 0;
    newTargets(t(source, ~switchEdge)) = 0;
    
    [~, ind] = sort(newTargets, 'descend');
    t(source, switchEdge) = ind(1:nnz(switchEdge));
end

h = graph(s,t);
Adj = adjacency(h);
end

function Adj = ...
    Graph_Generator_RenyiConnected(N_node,p_extralink)
% Erdos-Renyi network: each link (i,j) occurs with propability p
% Then ehanced to strong connectivity (random pick an out-neighbor if a node do not have out-neighbor in ER graph)
Adj = (rand(N_node) < p_extralink);

for i = 1:N_node
    Adj(i,i) = 0;
    for j = 1:i-1
        Adj(i,j) = Adj(j,i);
    end
end

for i = 1:N_node
    d_out = sum(Adj(i,:));
    if d_out == 0
        out_neighbor = randi(N_node)-1;
        if out_neighbor >= i
           out_neighbor = out_neighbor +1; 
        end
        Adj(i,out_neighbor) = 1;
        Adj(out_neighbor,i) = 1;
    end
end
end

function Adj = ...
    Graph_Generator_GEANT()
% GEANT graph (see cacheNetwork by Stratis)
Adj = zeros(22);
Adj(1,2) = 1;
Adj(1,22) = 1;
Adj(2,3) = 1;
Adj(2,4) = 1;
Adj(2,16) = 1;
Adj(3,4) = 1;
Adj(4,13) = 1;
Adj(4,5) = 1;
Adj(5,22) = 1;
Adj(5,19) = 1;
Adj(5,8) = 1;
Adj(6,7) = 1;
Adj(6,15) = 1;
Adj(7,8) = 1;
Adj(8,9) = 1;
Adj(9,10) = 1;
Adj(10,11) = 1;
Adj(10,13) = 1;
Adj(11,12) = 1;
Adj(12,13) = 1;
Adj(13,14) = 1;
Adj(13,16) = 1;
Adj(14,15) = 1;
Adj(15,16) = 1;
Adj(15,17) = 1;
Adj(15,22) = 1;
Adj(16,18) = 1;
Adj(17,18) = 1;
Adj(18,19) = 1;
Adj(19,20) = 1;
Adj(19,21) = 1;
Adj(20,21) = 1;
Adj(21,22) = 1;

Adj = Adj + Adj.'; % make symmetric
end

function Adj = ...
    Graph_Generator_LHC()
% DECO: LHC (Large Hadron Collider) topology is a prominent data-intensive computing network for high energy physics applications shown in Figure. 5(d)
Nodelist = {'S1','S2','S3','S4','S5','S6','S7','S8','NBR','UCSD','FNL','VND','UFL','WSC','MIT','PRD'};
G = graph([]);
G = addnode(G, Nodelist); 
G = addedge(G,'NBR','S6',1);
G = addedge(G,'NBR','FNL',1);
G = addedge(G,'NBR','S4',1);
G = addedge(G,'UCSD','S6',1);
G = addedge(G,'UCSD','FNL',1);
G = addedge(G,'UCSD','S8',1);
G = addedge(G,'FNL','S6',1);
G = addedge(G,'FNL','S8',1);
G = addedge(G,'FNL','S1',1);
G = addedge(G,'VND','S6',1);
G = addedge(G,'VND','S3',1);
G = addedge(G,'VND','UFL',1);
G = addedge(G,'UFL','S5',1);
G = addedge(G,'UFL','S7',1);
G = addedge(G,'UFL','S3',1);
G = addedge(G,'WSC','S6',1);
G = addedge(G,'WSC','S1',1);
G = addedge(G,'MIT','S6',1);
G = addedge(G,'MIT','S1',1);
G = addedge(G,'MIT','S5',1);
G = addedge(G,'PRD','S2',1);
G = addedge(G,'S1','S2',1);
G = addedge(G,'S1','S5',1);
G = addedge(G,'S2','S3',1);
G = addedge(G,'S2','S4',1);
G = addedge(G,'S4','S5',1);
G = addedge(G,'S4','S6',1);
G = addedge(G,'S5','S6',1);
G = addedge(G,'S5','S7',1);
G = addedge(G,'S6','S8',1);
G = addedge(G,'S7','S8',1);
Adj = full(adjacency(G));
end

function [Phi] = ...
    FtoPhi(Adj,N_app,Ta,Da,InputRate,f,g,eps)
% convert flow to phi var.
N_node = length(Adj);
N_stage = sum(Ta) + N_app;
G = digraph(Adj);
%N_edge = G.numedges * 2;
 [sOut,tOut] = findedge( G );
Edge = [sOut tOut]; % list of all edges (directed), each row is the [node_i,node_j]
N_edge = size(Edge,1);
Sa = Ta + 1; % number of stages for each app

% first calculate t using f,g
Mat_FtoT = kron(ones(1,N_node),eye(N_node*N_stage)); % input f is ordered i,j,s, need to sum all j,s
t = Mat_FtoT * f; % ti = sum_j fji + ri + input g
for node = 1:N_node
    for app = 1:N_app
        for k = 0:Ta(app)
            s_index = sum(Sa(1:app-1))+k +1;
            t_pos = (node -1)*N_stage + s_index;  % stage (a,0)
            g_pos = (node -1)*N_stage + s_index-1; % g for g(a,k-1) 
            if k == 0 % stage (a0), include r
                t(t_pos) = t(t_pos) + InputRate(node,app);
            else % stage (ak), cinlcude g
                t(t_pos) = t(t_pos) + g(g_pos);
            end
        end
    end
end

% then calculate phi
length_Phi = N_node * (1+N_node) * N_stage; % column vector odered as i ,j ,(d,m) and j can be 0
%length_Phi_plus = N_node * N_node * N_task; % column vector
Phi = zeros(length_Phi,1);
%Phi_plus = zeros(length_Phi_plus,1);
for app = 1:N_app
    da = Da(app);
    for k = 0: Ta(app)
        s_index = sum(Sa(1:app-1))+k +1;
        for node_i = 1:N_node
            t_pos = (node_i-1)*N_stage + s_index;
            if ~( k == Ta(app) && node_i == da) % if final stage destination, all phi 0

                if t(t_pos) < eps
                    % if t = 0, assign all phij to j=0 if k is not Ta, else pick the shortest path to destination
                    if k ~= Ta(app)
                        Phi_i0_pos = (node_i-1)* (N_node+1) * N_stage + s_index;
                        Phi(Phi_i0_pos) = 1;
                    else
                        P = shortestpath(G,node_i,da);
                        out_node = P(2);
                        Phi_pos = (node_i -1)*(N_node+1)*N_stage + (out_node)*N_stage + s_index;
                        Phi(Phi_pos) = 1;
                    end
                else % if t > 0, assign regularly
                    g_pos = (node_i-1)*N_stage + s_index;
                    Phi_i0_pos = (node_i-1)* (N_node+1) * N_stage + s_index;
                    Phi(Phi_i0_pos) = g(g_pos) / t(t_pos);
                    for node_j = 1:N_node
                        f_pos = (node_i -1)*N_node*N_stage + (node_j-1)*N_stage + s_index;
                        Phi_pos = (node_i -1)*(N_node+1)*N_stage + node_j*N_stage + s_index;
                        Phi(Phi_pos) = f(f_pos) / t(t_pos);
                    end
                end
            end
        end
    end
end
end


function [Is_Success,f_init,g_init] = ...
    Init_Generator2_MILP(Adj,N_app,Ta,Da,InputRate,PackSize,Delay_type,Delay_para,CompCost_type,CompCost_para, LinkCap,CompCap,Is_Initial_LCOR)
% Improved initial state generator, reduce the memory size
% formulate the LP/MILP wrt the link, not (node,node)

N_node = length(Adj);
N_stage = sum(Ta) + N_app;
G = digraph(Adj);
%N_edge = G.numedges * 2;
 [sOut,tOut] = findedge( G );
Edge = [sOut tOut]; % list of all edges (directed), each row is the [node_i,node_j]
N_edge = size(Edge,1);
Sa = Ta + 1; % number of stages for each app

%length_f_minus = N_edge * N_task;
%length_f_plus = N_edge * N_task;
%length_g = N_node * N_task;
length_f = N_edge * N_stage;
length_g = N_node * N_stage;


UB = sum(sum(InputRate))*1.5; % upper bound for all possible computation flow
computation_cap = UB; %(optional) a sharp upperbound additional to the linear computation cost

length_x_LP = length_f + length_g; 
% dimension of x vetcor, arranged with f_e(a,k) forall stage, g_i(a,k) forall stage

length_constraints = N_node * N_stage + N_edge + N_node; 
% flow conserv for f at i. And LinkCap for all (i,j), CompCap for all i.

if ~Is_Initial_LCOR
    % if not requaired Initial LCOR, just random objective vector
    %f_LP = ones(1,length_x_LP);
    f_LP = rand(1,length_x_LP);
else
    % if required initial local, the objecteive is set to the following:
    % entries for f stage 0 are 1e3 (sufficiently large), to ensure the data flow is minimized;
    % entries for f other stage are given by the marginal cost of link ij at Fij=0, in order to generate a shortest path for result flow
    % entries for g are 0, as compute flow are not to be minimized
    f_LP = zeros(1,length_x_LP);
    for node_i = 1:N_node
        for node_j = find(Adj(node_i,:)) % specify the marginal at F=0
            if strcmp(Delay_type,'queue') % if queueing delay D = F/(C-F), D'(0) = C/(C-F)^2 = 1/C
                Margin_Dij = 1 / Delay_para(node_i,node_j);
            elseif strcmp(Delay_type,'linear') % if linear cost
                Margin_Dij = Delay_para(node_i,node_j);
            else
                error('ERROR: unknown delay type\n');
            end
            for app = 1:N_app
                for k = 0:Ta(app) % all stages
                    s_index = sum(Sa(1:app-1))+k +1; % stages are counted 0-Ta for each app 
                    [tf, edge_index]=ismember([node_i node_j],Edge,'rows');
                    f_pos = (edge_index -1)*N_stage + s_index;  % f is pderdered edge, stage
                    if k== 0
                        f_LP(f_pos) = 1e3; % Shortest path with marginal at F=0
                    else
                        f_LP(f_pos) = Margin_Dij * PackSize(s_index);
                    end
                end
            end
        end
    end
end

A_LP_eq = zeros(N_node * N_stage,length_x_LP);
b_LP_eq = zeros(N_node * N_stage,1);
A_LP_neq = zeros( N_edge + N_node,length_x_LP);
b_LP_neq = zeros( N_edge + N_node,1);

for node_i = 1:N_node
    % first the flow conservarions
    % note: spacial case at k = 0, k = Ta and k = da are not included, they're included in the UB
    for app = 1:N_app
        for k = 0:Ta(app)
            s_index = sum(Sa(1:app-1))+k +1;
            % position in constraints
            index_conserv = (node_i-1)*N_stage + s_index; % constraints are orderd node, stage

            % set the conserv for f
            % \sum_{j}fij + gi(k) = sum_j fji + gi(k-1) + r_i
            if k == 0 % r for k = 0
                b_LP_eq(index_conserv) = -1* InputRate(node_i,app);
            else
                b_LP_eq(index_conserv) = 0;
            end

            index_g = length_f + (node_i-1)*N_stage + s_index;
            if k ~= 0
                A_LP_eq(index_conserv,index_g-1) = 1; % results from stage k-1
            end
            if k ~= Ta(app)
                A_LP_eq(index_conserv,index_g) = -1; % output to CPU
            end

            for node_j = 1:N_node
                if ~( k == Ta(app) && node_i == Da(app))
                    if Adj(node_j,node_i) > 0 % in-going links (j,i)
                        [tf, edge_index]=ismember([node_j node_i],Edge,'rows');
                        index_f = (edge_index -1)*N_stage + s_index;
                        A_LP_eq(index_conserv,index_f) = 1;
                    end
                end

                if Adj(node_i,node_j) > 0 % out-going links (i,j)
                    [tf, edge_index]=ismember([node_i node_j],Edge,'rows');
                    index_f = (edge_index -1)*N_stage + s_index;
                    A_LP_eq(index_conserv,index_f) = -1;
                end

            end

        end
    end

    % then the link capacity constraint for all (i,j), sum_s fijs <= LinkCap
    for node_j = find(Adj(node_i,:)) % for all (i,j) in E
        [tf, edge_index]=ismember([node_i node_j],Edge,'rows');
        index_LinkCap =  edge_index; % constraint position
        for s_index = 1:N_stage
            index_f =  (edge_index -1)*N_stage + s_index;
            A_LP_neq(index_LinkCap,index_f) = 1;
        end
        b_LP_neq(index_LinkCap) = LinkCap(node_i,node_j);
    end

    % then the computation capacity for all nodes, sum_t git <= CompCap
    index_CompCap = N_edge + node_i;
    for s_index = 1:N_stage
        index_g = length_f + (node_i -1)*N_stage + s_index;
        A_LP_neq(index_CompCap,index_g) = 1;
    end
    b_LP_neq(index_CompCap) = CompCap(node_i);
end

x_lb_LP = zeros(1,length_x_LP);
x_ub_LP = UB * ones(1,length_x_LP);
% restrictions on f and g for special cases
% 1. if k = 0, g flows from k-1 does not exist (already in A)
% 2. if k = Ta, g to k should be 0
% for node_i = 1:N_node
%     % for node_j = find(Adj(node_i,:)) % for all (i,j) in E
%     %    [tf, edge_index]=ismember([node_i node_j],Edge,'rows');
%     for app = 1:N_app
%         k = Ta(app);
%         s_index = sum(Sa(1:app-1))+k +1;
%         index_g = length_f + (node_i-1)*N_stage + s_index;
%         x_ub_LP(index_g) = 0;
%     end
%     %end
% end
% 3. if k = Ta and i = da, out-going f should be 0
% for app = 1:N_app
%     da = Da(app);
%     k = Ta(app);
%     s_index = sum(Sa(1:app-1))+k +1;
%     for node_j = find(Adj(da,:)) % for all (i,j) in E
%         [tf, edge_index]=ismember([da node_j],Edge,'rows');
%         index_f =  (edge_index -1)*N_stage + s_index;
%         x_ub_LP(index_f) = 0;
%     end
% end

options =  optimoptions(@linprog,'Display','off');
fprintf('Init_Generator2: begin solving LP...\n');
[x_opt_LP, cost_opt_LP,exitflag] = linprog(f_LP,A_LP_neq,b_LP_neq,A_LP_eq,b_LP_eq,x_lb_LP,x_ub_LP);
%[x_opt_LP, cost_opt_LP,exitflag] = linprog(f_LP,A_LP,b_LP,[],[],x_lb_LP,x_ub_LP,options);
fprintf('Init_Generator: LP finished\n');

if exitflag <= 0 % No feasible solution or unbounded
    Is_Success = 0;
    f_init_edge = zeros(length_f,1);
    g_init_edge = zeros(length_g,1);
else
    Is_Success = 1;
    x_opt_LP = reshape(x_opt_LP,[],1);
    f_init_edge = x_opt_LP(1:length_f);
    g_init_edge = x_opt_LP(length_f + 1: length_f + length_g);
end

% finally change edge-based variable to node pair-based
length_f_output = N_node*N_node*N_stage;
length_g_output = N_node * N_stage;
f_init = zeros(length_f_output,1);
%g_comp_init = zeros(length_g_output,1);
for s_index = 1:N_stage
    for edge_index = 1:N_edge
        node_i = Edge(edge_index,1);
        node_j = Edge(edge_index,2);
        pos_f_edge = (edge_index -1)*N_stage + s_index;
        pos_f_init = (node_i -1)*N_node*N_stage + (node_j-1)*N_stage + s_index;
        f_init(pos_f_init) = f_init_edge(pos_f_edge);
    end
end
g_init = g_init_edge;

end

function [Fij,Gi,t,Is_Loopfree] = Update_Flow(Adj,N_app,Ta,Da,InputRate,PackSize,CompWeight,Phi,eps)
% calculate flow Fij, Gi and ti from given phi
% Fij: matrix (N x N)
% Gi: vector (N x 1)
% t: column vector (N*N_stage, 1)
% CompWeight: matrix (N_node x N_stage)
% Adj: adj matrix; M: computation types; Phi is given in column vector

N_node = length(Adj);
N_stage = sum(Ta) + N_app;
G = digraph(Adj);
%N_edge = G.numedges * 2;
 [sOut,tOut] = findedge( G );
Edge = [sOut tOut]; % list of all edges (directed), each row is the [node_i,node_j]
N_edge = size(Edge,1);
Sa = Ta + 1; % number of stages for each app

MAX_PATH_LEN = N_node;
Is_Loopfree = 1; % if the input phi is loop-free

length_Phi = N_node * (1+N_node) * N_stage; % column vector odered as i ,j ,(d,m) and j can be 0
length_f = N_node * N_node * N_stage; % column vector ordered as i, j, (d,m)
length_g = N_node * N_stage; % column vector ordered as i, (d,m)
length_t = N_node * N_stage; % column vector

% Method: iterate calculate ti = \sum j tj * phi_j + r + g, start from data source, run for at most N
% times, recursive for stages from Ta to 0
t = zeros(length_t,1);
f = zeros(length_f,1);
g = zeros(length_g,1);
Fij = zeros(N_node,N_node); % Fij for data, matrix
Gi = zeros(N_node,1); % Gi

for app = 1:N_app
    da = Da(app);
   
    % NEW METHOD: solve the set of linear equations: 
    % t_i - \sum_j t_j phi_ji = r + g_input
    % solve Ax=b,   

    for k = 0: Ta(app)  % recursive calculate from final stage to
        s_index = sum(Sa(1:app-1))+k +1;

        % first test loop-free
        DAG_stage = zeros(N_node,N_node); % Adj matrix for directed graph of current stage
        for node_i = 1:N_node
            for node_j = 1:N_node
                Phi_pos = (node_i -1)* (N_node+1) * N_stage + node_j * N_stage + s_index;
                if Phi(Phi_pos) > eps % mark the link with phi > 0
                    if Adj(node_i,node_j) == 0
                        error('ERROR: Positive phi+ on non-link (%d,%d)!\n',node_i,node_j)
                    end
                    DAG_stage(node_i,node_j) = Phi(Phi_pos);
                end
            end
        end
        G_s = digraph(DAG_stage);
        Is_Loopfree = 1 - hascycles(G_s);
        if Is_Loopfree == 0
            pause();
        end

        % extract phi into matrix
        Phi_s = zeros(N_node,N_node); % phi- for current task, matrix, exclude phi_i0
        Phi_0_s = zeros(N_node,1); % phi_i0 for current task, column vector
        for node_i = 1:N_node
            phi_i0_pos = (node_i-1) * (1+N_node)*N_stage + s_index;
            Phi_0_s(node_i) = Phi(phi_i0_pos);
            for node_j = 1:N_node
                phi_pos = (node_i-1) * (1+N_node) * N_stage + node_j * N_stage + s_index; % Note to avoid phi_i0
                Phi_s(node_i,node_j) = Phi(phi_pos);
            end
        end
        
        % calculate r+g_input
        inputVec = zeros(N_node,1);
        if k == 0
            inputVec = InputRate(:,app);
        else
            for node_i = 1:N_node
                inputVec(node_i) = g((node_i-1) * N_stage + s_index -1);
            end
        end
        % varialbe length
        len_x_lin = N_node;
        len_b_lin = N_node; % number of constraint
        % matrix A,b
        b_lin = inputVec;
        A_lin = zeros(len_b_lin,len_x_lin);
        for node_i = 1:N_node
            A_lin(node_i,node_i) = 1;
            for node_j = find(Adj(:,node_i)) % j,i in E
                A_lin(node_i,node_j) = -1 * Phi_s(node_j,node_i);
            end
        end
        t_s = linsolve(A_lin,b_lin);

        % calculate g
        g_s = t_s .* Phi_0_s; % gidm, column vector

        % assign entries
        for node_i = 1:N_node
            g_pos = (node_i-1)*N_stage + s_index;
            g(g_pos) = g_s(node_i);
            t_pos = (node_i-1)*N_stage + s_index;
            t(t_pos) = t_s(node_i);
            Gi(node_i) = Gi(node_i) + g_s(node_i) * CompWeight(node_i,s_index);
            for node_j = 1:N_node
                f_pos = (node_i-1)*N_node*N_stage + (node_j-1)*N_stage + s_index;
                f(f_pos) = t_s(node_i) * Phi_s(node_i,node_j);
                Fij(node_i,node_j) = Fij(node_i,node_j) + f(f_pos) * PackSize(s_index);    
            end
        end
    end
end

end

function [Delay,Delay_deri,CompCost,CompCost_deri] = Update_Cost(Delay_type,CompCost_type,Delay_para,CompCost_para,Fij,Gi,eps)
% delay calculation. Fij: N_node * N_node;
% Delay: N * N; Delay_deri: N * N
% Computation cost calculation, Comp_Flow: (N_node x1) ; CompCost: N * 1; CompCost_deri: N * 1
% Queueing Delay: D = F/(C-F)
N_node = length(Fij);
Delay = zeros(size(Fij));
Delay_deri = zeros(size(Fij));
if strcmp(Delay_type,'queue')
    Cap =  Delay_para;
    for node_i = 1:N_node
        for node_j = 1:N_node
            c = Cap(node_i,node_j);
            f = Fij(node_i,node_j);
            if c < eps % if there is no link ij
                if f < eps % if no flow
                    Delay(node_i,node_j) = 0;
                    Delay_deri(node_i,node_j) = 0;
                else
                    Delay(node_i,node_j) = Inf;
                    Delay_deri(node_i,node_j) = Inf;
                    error('ERROR: Positive flow on non-link');
                end
            else % if there is link ij
                if c-f < eps % if exceed capacity
                    Delay(node_i,node_j) = Inf;
                    Delay_deri(node_i,node_j) = Inf;
                    error('ERROR: Exceed Link Capacity');
                else
                    Delay(node_i,node_j) = f/(c-f);
                    Delay_deri(node_i,node_j) = c/(c-f)^2;
                end
            end
        end
    end
elseif strcmp(Delay_type,'linear')
    Delay = Delay_para .* Fij;
    Delay_deri = Delay_para;
end

if strcmp(CompCost_type,'sum_queue')
    Cap =  CompCost_para;
    if min(Cap - Gi) <= 0
        %CompCost = Inf * Comp_Flow_Sum;
        %CompCost_deri = Inf * CompFlow;
        error('ERROR: Exceed Computation Capacity');
    else
        CompCost = Gi ./ (Cap - Gi);
        CompCost_deri = Cap ./ (Cap - Gi).^2;
    end
elseif strcmp(CompCost_type,'sum_linear')
    CompCost = CompCost_para .* Gi;
    CompCost_deri = CompCost_para;
    %CompCost_deri = kron(CompCost_deri,ones(1,size(CompFlow,2)));
end

end

function Is_Valid = Check_Phi_Valid(Adj,N_app,Ta,Da,InputRate,Phi,PackSize,CompWeight,LinkCap, CompCap, eps)
% check if given phi is valid:
% 1: are flow conervation satisfied
% 2: are + / - flow loop-free
% 3: are capacity constraints satisfied
Is_ShowDetail = 1;
Is_Valid = 1;
N_node = length(Adj);
N_stage = sum(Ta) + N_app;
G = digraph(Adj);
%N_edge = G.numedges * 2;
 [sOut,tOut] = findedge( G );
Edge = [sOut tOut]; % list of all edges (directed), each row is the [node_i,node_j]
N_edge = size(Edge,1);
Sa = Ta + 1; % number of stages for each app

% Step 1: all phi sums up to 1, 
% except final stage g = 0; final stage destination sum up to 0  
for app = 1:N_app
    da = Da(app);
    for k = 0:Ta(app)
        s_index = sum(Sa(1:app-1))+k +1;
        for node_i = 1:N_node
            sum_phi = 0;
            for node_j = [0 find(Adj(node_i,:))]
                Phi_pos = (node_i-1)*(1+N_node)*N_stage + node_j*N_stage + s_index;
                sum_phi = sum_phi + Phi(Phi_pos);
                if k == Ta(app) && node_j == 0 && Phi(Phi_pos) > eps
                    if Is_ShowDetail
                        fprintf('Validation Fail: Positive comp flow for final stage.\n');
                    end
                    Is_Valid = 0;
                    return
                end
            end
            if k == Ta(app) && node_i == da
                if abs(sum_phi) >= eps % if is detinaton and sum of phi+ is not 0
                    if Is_ShowDetail
                        fprintf('Validation Fail: Sum of phi is not 0 at final stage destination.\n');
                    end
                    Is_Valid = 0;
                    return
                end
            else
                if abs(sum_phi - 1) >= eps % if is not dest and sum of phi+ is not 1
                    if Is_ShowDetail
                        fprintf('Validation Fail: Sum of phi is not 1');
                    end
                    Is_Valid = 0;
                    return
                end
            end
           
        end
    end
end

% Step 2: compute the flow, while check loop-free
[Fij,Gi,t,Is_Loopfree] = ...
    Update_Flow(Adj,N_app,Ta,Da,InputRate,PackSize,CompWeight,Phi,eps);
if Is_Loopfree == 0
    if Is_ShowDetail
        fprintf('Validation Fail: Contain loop.\n');
    end
    Is_Valid = 0;
    return
end

% Step 3: check capacity
if min(LinkCap - Fij) <= -eps
    if Is_ShowDetail
        fprintf('Validation Fail: Exceed link capacity.\n');
    end
    Is_Valid = 0;
    return
end
if min(CompCap - Gi) <= -eps
    if Is_ShowDetail
        fprintf('Validation Fail: Exceed computation capacity.\n');
    end
    Is_Valid = 0;
    return
end

end

function [TotalCost_Opt_Offline, f_minus_Opt_offline, f_plus_Opt_offline,g_Opt_offline] = ...
    Offline_Optimization(Adj,Task,InputRate,a_equivalent,Delay_type,CompCost_type,Delay_para,CompCost_para,...
    Link_Cap,Comp_Cap,f_minus_init,f_plus_init,g_init,Capacity_Saturate_Factor,eps)
% Flow model centralized solver for global optimal with given start point, with b=0.
% a_equivalent is not meant for every m but for every i,(d,m), to incorporate lower bound calculation
% a_equivalent: column vector (N*N_task x 1)
% TotalCost_Opt_Offline: scalar
% f_minus_Opt_offline, f_plus_Opt_offline: column vector (N*N*N_task x 1)
% g_Opt_offline: column vector (N*N_task x 1)
% method: fmincon with specified gradient and Hessian
% x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
N_node = length(Adj);
N_task = size(Task,1);
length_x = N_node * N_node * N_task *2 + N_node * N_task; % length of variable, ordered as f-ijdm, f+ijdm, gidm

% construct constraints
length_Aeq = 2 * N_node * N_task; % flow conservation for f-, f+ at each node for each task, ordered as f-idm, f+idm
length_Aleq = N_node * N_node + N_node; % link flow upperbound and node computation upperbound (capacity * Capacity_Saturate_Factor)
Aeq_fmin = zeros(length_Aeq,length_x);
beq_fmin = zeros(length_Aeq,1);
Aleq_fmin = zeros(length_Aleq,length_x);
bleq_fmin = zeros(length_Aleq,1);
lb_fmin = zeros(length_x,1);
ub_fmin = zeros(length_x,1); % note: ub to be determind 0 or Inf by topology

% first flow conservation
for t_index = 1:N_task
    t_dest = Task(t_index,1);
    for node_i = 1:N_node
        
        % data flow conservation: (sumj f-jidm) - (sumj f-ijdm + gidm)  <= -ridm
        Aeq_minus_pos = (node_i -1)*N_task + t_index;
        ridm = InputRate(node_i,t_index);
        gidm_pos = 2*N_node*N_node*N_task + (node_i -1)*N_task + t_index; % position of g in x
        for node_j = find(Adj(node_i,:)) % fijdm-
            fijdm_minus_pos = (node_i-1)*N_node*N_task + (node_j-1)*N_task + t_index; % pos of fijdm- in x
            Aeq_fmin(Aeq_minus_pos,fijdm_minus_pos) = -1;
        end
        Aeq_fmin(Aeq_minus_pos,gidm_pos) = -1;
        for node_j = find(Adj(:,node_i)) % fjidm-
            fjidm_minus_pos = (node_j-1)*N_node*N_task + (node_i-1)*N_task + t_index;
            Aeq_fmin(Aeq_minus_pos,fjidm_minus_pos) = 1;
        end
        beq_fmin(Aeq_minus_pos) = -ridm;
        
        % result flow conservation: -(sumj f+ijdm) <= 0 if is dest, and (a*gidm + sumj f+jidm) - (sumj f+ijdm) <= 0 if not dest
        Aeq_plus_pos = N_node*N_task + (node_i-1)*N_task + t_index;
        a_pos = (node_i-1)*N_task + t_index;
        for node_j = find(Adj(node_i,:)) % fijdm+
            fijdm_plus_pos = N_node*N_node*N_task + (node_i-1)*N_node*N_task + (node_j-1)*N_task + t_index;
            Aeq_fmin(Aeq_plus_pos,fijdm_plus_pos) = -1;
        end
        if node_i ~= t_dest
            Aeq_fmin(Aeq_plus_pos,gidm_pos) = a_equivalent(a_pos);
            for node_j = find(Adj(:,node_i)) % fjidm+
                fjidm_plus_pos = N_node*N_node*N_task + (node_j-1)*N_node*N_task + (node_i-1)*N_task + t_index;
                Aeq_fmin(Aeq_plus_pos,fjidm_plus_pos) = 1;
            end
        end
    end
end

% then link flow and computation capacity
for node_i = 1:N_node
    for node_j = 1:N_node % link capacities: sum_t f-ijt + f+ijt <= Cij
        Aleq_link_cap_pos = (node_i-1)*N_node + node_j;
        for t_index = 1:N_task
            fijdm_minus_pos = (node_i-1)*N_node*N_task + (node_j -1)*N_task + t_index;
            fijdm_plus_pos = N_node*N_node*N_task + (node_i-1)*N_node*N_task + (node_j -1)*N_task + t_index;
            Aleq_fmin(Aleq_link_cap_pos,fijdm_minus_pos) = 1;
            Aleq_fmin(Aleq_link_cap_pos,fijdm_plus_pos) = 1;
        end
        bleq_fmin(Aleq_link_cap_pos) = Link_Cap(node_i,node_j) * Capacity_Saturate_Factor;
    end
    Aleq_comp_cap_pos = N_node*N_node + node_i;
    for t_index = 1:N_task % computation capacity: sum_t git <= Ci
        g_idm_pos = N_node*N_node*N_task*2 + (node_i-1)*N_task + t_index;
        Aleq_fmin(Aleq_comp_cap_pos,g_idm_pos) = 1;
    end
    bleq_fmin(Aleq_comp_cap_pos) = Comp_Cap(node_i) * Capacity_Saturate_Factor;
end

% then upper bound of f-,f+,g according to topology
MAX_FLOW = sum(sum(InputRate)) * max(1,max(a_equivalent));
for t_index = 1:N_task
    for node_i = 1:N_node
        for node_j = 1:N_node
            f_minus_pos = (node_i-1)*N_node*N_task + (node_j-1)*N_task + t_index;
            f_plus_pos = N_node*N_node*N_task + (node_i-1)*N_node*N_task + (node_j-1)*N_task + t_index;
            if Adj(node_i,node_j) == 1
                ub_fmin(f_minus_pos) = MAX_FLOW;
                ub_fmin(f_plus_pos) = MAX_FLOW;
            else
                ub_fmin(f_minus_pos) = 0;
                ub_fmin(f_plus_pos) = 0;
            end
        end
        g_pos = 2*N_node*N_node*N_task + (node_i-1)*N_task + t_index;
        ub_fmin(g_pos) = MAX_FLOW;
    end
end

% construct objective functions wrt different cost type, with output gradient
length_Flow = N_node*N_node + N_node; % Flow: [Fij Gi]
x_to_Flow = zeros(length_Flow,length_x); % calculate [Fij Gi] from x: Flow = x_to_Flow * x
Flow_Para = zeros(length_Flow,1); % vector of parameter [Cij Ci]
for node_i = 1:N_node
    for node_j = find(Adj(node_i,:))
        Flow_pos = (node_i-1)*N_node + node_j;
        for t_index = 1:N_task
            f_minus_pos = (node_i-1)*N_node*N_task + (node_j-1)*N_task + t_index;
            f_plus_pos = N_node*N_node*N_task + (node_i-1)*N_node*N_task + (node_j-1)*N_task + t_index;
            x_to_Flow(Flow_pos,f_minus_pos) = 1;
            x_to_Flow(Flow_pos,f_plus_pos) = 1;
        end
        Flow_Para(Flow_pos) = Delay_para(node_i,node_j);
    end
    Flow_comp_pos = N_node*N_node + node_i;
    for t_index = 1:N_task
        g_pos = N_node*N_node*N_task*2 + (node_i-1)*N_task + t_index;
        x_to_Flow(Flow_comp_pos,g_pos) = 1;
    end
    Flow_Para(Flow_comp_pos) = CompCost_para(node_i);
end

Valid_Link_Pos_List = find(Flow_Para > eps); % mark the positions of valid link (link capa/computation cap is positive)
x_to_Flow_Valid = x_to_Flow(Valid_Link_Pos_List,:);
Flow_Para_Valid = Flow_Para(Valid_Link_Pos_List);

    function [TotalCost,Gradient] = Cost_Queue(x)
        % total cost sum_ij Dij(Fij) + sumi Ccompi(Gi),
        % Fij = sumt f-ijt + f+ijt, Gi = sumt git, Dij(y) = y/(Cij - y) , Ccompi(y) = c/(Ci-y)
        TotalCost = sum( (x_to_Flow_Valid * x)./(Flow_Para_Valid- x_to_Flow_Valid * x) );
        % Gradient : D'ij = C/(C-F)^2, pD/pfij- = D'ij
        if nargout > 1
            Flow = x_to_Flow * x;
            Delay_grad = (Flow_Para + eps) ./ ((Flow_Para -  Flow).^2 + eps); % add eps to avoid Inf and NaN
            Delay_grad_augmented = [Delay_grad(1:N_node*N_node) ; Delay_grad]; % [D'ij; D'ij; G'i]
            Gradient = kron( Delay_grad_augmented, ones(N_task,1) ); % [D'ij x t; D'ij x t; G'i x t ]
        end
    end
    function [TotalCost,Gradient] = Cost_Linear(x)
        % total cost sum_ij Dij(Fij) + sumi Ccompi(Gi),
        % Fij = sumt f-ijt + f+ijt, Gi = sumt git, Dij(y) = y/(Cij - y) , Ccompi(y) = c/(Ci-y)
        TotalCost = sum((x_to_Flow_Valid * x) .* Flow_Para_Valid);
        % Gradient
        if nargout > 1
            Gradient = kron([Flow_Para(1:N_node*N_node);Flow_Para],ones(N_task,1)); % [Cij;Cij;Ccompi]
        end
    end

% solve fmincon
x0 = [f_minus_init;f_plus_init;g_init];
[TotalCost_init,Cost_gradient_init] = Cost_Queue(x0);
options =  optimoptions(@fmincon,'Display','off','MaxFunctionEvaluations',3e4, 'MaxIterations',1e4,'SpecifyObjectiveGradient',true);
options =  optimoptions(@fmincon,'Display','iter-detailed','MaxFunctionEvaluations',1e3, 'MaxIterations',1e2,'SpecifyObjectiveGradient',true);
if strcmp(Delay_type,'queue') && strcmp(CompCost_type,'sum_queue') % if use queueing delay for both computation and link
    %[x_opt,f_val,exitflag] = fmincon(@Cost_Queue,x0,Aleq_fmin,bleq_fmin,Aeq_fmin,beq_fmin,lb_fmin,ub_fmin,[],options);
    [x_opt,f_val,exitflag] = fmincon(@Cost_Queue,x0,[Aleq_fmin; -eye(length_x)],[bleq_fmin; zeros(length_x,1)],Aeq_fmin,beq_fmin,[],[],[],options);
elseif strcmp(Delay_type,'linear') && strcmp(CompCost_type,'sum_linear') % if use linear cost for both computation and link
    [x_opt,f_val,exitflag] = fmincon(@Cost_Linear,x0,Aleq_fmin,bleq_fmin,Aeq_fmin,beq_fmin,lb_fmin,ub_fmin,[],options);
else
    error('ERROR: wrong cost type!\n');
end

if exitflag < 0
    error('ERROR: fail to solve.\n');
elseif exitflag == 0
    error('ERROR: exceed max iteration.\n');
end
TotalCost_Opt_Offline = f_val;
f_minus_Opt_offline = x_opt(1 : N_node*N_node*N_task);
f_plus_Opt_offline = x_opt(N_node*N_node*N_task+1 : 2*N_node*N_node*N_task);
g_Opt_offline = x_opt(2*N_node*N_node*N_task +1: end);
%f_minus_Opt_offline_mat = reshape(f_minus_Opt_offline,N_node,N_node)';
%f_plus_Opt_offline_mat = reshape(f_plus_Opt_offline,N_node,N_node)';
end

function [Lifetime_data,Lifetime_result] = Average_Liftime(LinkFlow_data,LinkFlow_result,Delay,Task,InputRate,a_m)
% calculate average package lifetime for a data/result package, given D_ij is the package delay on (i,j)
% L_data = (\sum Dij*F-ij)/(\sum ridm); L_result = (\sum Dij*F+ij)/(\sum am*ridm) 
% Note: b is not considered, since b is extra flow that should be minimized
TotalCost_data = sum(sum(LinkFlow_data .* Delay));
TotalCost_result = sum(sum(LinkFlow_result .* Delay));
TotalRate_data = sum(sum(InputRate));
TotalRate_result = 0;
for node_i = 1:size(InputRate,1)
    for t_index = 1:size(InputRate,2)
        t_m = Task(t_index,2);
        TotalRate_result = TotalRate_result + InputRate(node_i,t_index) * a_m(t_m);
    end
end
Lifetime_data = TotalCost_data / TotalRate_data;
Lifetime_result = TotalCost_result / TotalRate_result;
end



