% Flow model simulator for the joint routing/compuation with general convex cost.
% Created by Jinkun, Feb 20 2022.
clear all
close all
clc

Is_Save = 0; % 1: generated parameter will be saved
Is_Load = 1; % 1: will abandon generated parameter and load saved file
Is_debugTest = 0; % 1: debug test for first iteration
Network_Type = 'random_thread';
eps = 1e-5; % threshold for justify 0
%Network_Type = 'example_1'; % example in figure 3
%% Graph Generating
N_node = 10; % number of nodes

if strcmp(Network_Type,'random_thread')
    p_extralink = 0.2; % graph is a linear network with extra random link
    Adj = Graph_Generator_RandomThread(N_node,p_extralink);
elseif strcmp(Network_Type,'example_1')
    N_node = 4;
    Adj = [0 1 0 1; 1 0 1 1; 0 1 0 1; 1 1 1 0];
end

%% Task and rates
Task_Type = 'random'; % random tasks
N_task = 7;
M = 3; % number of computation type

% Task: matrix (N_task x 2 ), each row is (d,m)
% InputRate: matrix (N_node x N_task), each row is input rate for all task of one node
if strcmp(Task_Type,'random')
    inpute_nodes_number_max = 5; % the maximum input node number of each task
    rate_min = 1;
    rate_max = 2; % max input rate of task (uniform)
    [Task,InputRate] = Task_Generator_random(Adj,N_task,M,rate_min,rate_max,inpute_nodes_number_max);
end
if strcmp(Network_Type,'example_1')
    N_task = 2;
    M = 1;
    a_m = 0.5;
    b_overhead = 0.2;
    Task = [4 1];
    InputRate = [1 0 0 0]';
end

%% phi variables and convert functions
% calculate flow Fij, Gim and Ti+- from routing variable phi, and reverse
% for phi -> f, see function 'Update_Flow' in appendix
% for f -> phi, see function 'FtoPhi' in appendix

%% Delay and Derivative
Delay_type = 'queue'; % queueing delay F/(C-F)
%Delay_type = 'linear'; % linear link cost

CompCost_type = 'sum_queue'; % sum of all m and apply queueing delay
%CompCost_type = 'sum_linear'; % sum of all m and linear

if strcmp(Network_Type,'random_thread')
    Cap_Max = 10;
    Cap_Min = 5;
    Delay_para = (Cap_Min + (Cap_Max-Cap_Min) * rand(N_node)) .* Adj; % parameter for delay, e.g. link capacity. Matrix (N x N)
    CompCost_para = Cap_Min + (Cap_Max-Cap_Min) * rand(N_node,1); % parameter for computation cost, e.g. CPU speed. Column vecror (N x 1)
elseif strcmp(Network_Type,'example_1')
    Delay_para = 2*Adj; % all link capa are 2
    CompCost_para = [1 1 1 2]';
end
% see function 'Update_Cost' in appendix


%% Save and load topology, parameters
SaveTopoFileName = 'RouteComputeNetwork_SaveTopoPara';
if Is_Save
    save(SaveTopoFileName);
end
if Is_Load
    load(SaveTopoFileName);
end
figure(1)
G = graph(Adj);
plot(G)

%% data-result convert rate
a_m = rand(1,M);
a_m = 0.5*ones(1,M);
b_overhead = 0; % b in paper

%% Generate initial state
Initial_Type = 'MILP'; % find a feasibe point using MILP with some arbitrary objective

if strcmp(Network_Type,'example_1')
    Phi_minus_init = [0 0 0 0 1 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0];
    Phi_plus_init = [0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0];
    
else
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
    if strcmp(Initial_Type,'MILP')
        Capacity_Saturate_Factor = 0.9; % the maximum allowed fraction of capacity
        [Is_Success,f_minus_init,f_plus_init,g_comp_init] = Init_Generator_MILP(Adj,Task,InputRate,a_m,b_overhead,LinkCap*Capacity_Saturate_Factor,CompCap*Capacity_Saturate_Factor);
        if ~Is_Success
            error('ERROR: No Feasible Initial State!  Will abort.\n');
        end
        [Phi_minus_init,Phi_plus_init] = FtoPhi(Adj,Task,InputRate,a_m,b_overhead,f_minus_init,f_plus_init,g_comp_init,eps);
    end
    if strcmp(Network_Type,'example_1')
        Phi_minus_init = [0 0 0 0 1 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0];
        Phi_plus_init = [0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0];
    end
end
%% Algorithm Parameters
T_MAX = 100; % max iteration number;
EPS_STOP_ratio = 1e-3; % stop if total cost does not decrease more than this ratio

Is_Use_GP = 1;      % if using method: Gradient Projection (scaled matrix = I)
Is_Use_SGP = 0;     % if using method: Scaled Gradient Projection
Is_Use_RA = 0;      % if using method: random allocation
Is_Use_1Hop = 0;    % if using method: optimal 1-hop
Is_Use_Local = 0;   % if using method: local computation (with necessary offloading)
Is_Use_BP = 0;      % if using method: Backpressure
Is_Use_LB = 0;      % if using method: Lower bound generalized with centralized flow-model optimization with am replaced by am+b/ti-, not accurate
Is_Use_Offline = 1;     % if using method: contralized flow-model optimization that ignore b, only performed at initailziation, not accurate

if Is_Use_GP
    %StepSize_GP = 0.1./ sqrt(1:T_MAX); % use universal stepsize alpha = c/sqrt(t);
    StepSize_GP = 0.1* ones(1,T_MAX); % use universal stepsize alpha = c/sqrt(t);
    TotalCost_GP = NaN * ones(1,T_MAX);
    Phi_minus_GP_current = Phi_minus_init;
    Phi_plus_GP_current = Phi_plus_init;
end
if Is_Use_SGP
    TotalCost_SGP = NaN * ones(1,T_MAX);
    Phi_minus_SGP_current = Phi_minus_init;
    Phi_plus_SGP_current = Phi_plus_init;
end
if Is_Use_LB
    TotalCost_LB = NaN * ones(1,T_MAX);
    Phi_minus_LB_current = Phi_minus_init;
    Phi_plus_LB_current = Phi_plus_init;
end
if Is_Use_Offline
    % generate equivalent a_idm
    a_equivalent = zeros(N_node*N_task,1);
    for t_index = 1:N_task
        t_m = Task(t_index,2);
        for node_i = 1:N_node
            a_equivalent((node_i-1)*N_task + t_index) = a_m(t_m);
        end
    end
    [TotalCost_Opt_Offline, f_minus_Opt_offline, f_plus_Opt_offline,g_Opt_offline] = ...
        Offline_Optimization(Adj,Task,InputRate,a_equivalent,Delay_type,CompCost_type,Delay_para,CompCost_para,...
        LinkCap,CompCap,f_minus_init,f_plus_init,g_comp_init,Capacity_Saturate_Factor,eps);
    TotalCost_Offline = TotalCost_Opt_Offline * ones(1,T_MAX);
end

%% run

for iter_t = 1:T_MAX
    if mod(iter_t, 1 ) == 0
        fprintf('Iteration no. %d\n',iter_t);
    end
    if Is_Use_GP % if using gradient projection algorithm
        % update and save current statues
        [LinkFlow_GP_current,CompFlow_GP_current,t_minus_GP_current,t_plus_GP_current] ...
            = Update_Flow(Adj,M,Task,a_m,b_overhead,InputRate,Phi_minus_GP_current,Phi_plus_GP_current,eps);
        
        [Delay_GP_current,Delay_deri_GP_current,CompCost_GP_current,CompCost_deri_GP_current] ...
            = Update_Cost(Delay_type,CompCost_type,Delay_para,CompCost_para,LinkFlow_GP_current,CompFlow_GP_current,eps);
        
        TotalCost_GP(iter_t) = sum(sum(Delay_GP_current)) + sum(sum(CompCost_GP_current));
        %         if iter_t >= 2
        %             if  abs(TotalCost_GP(iter_t) - TotalCost_GP(iter_t -1))/TotalCost_GP(iter_t -1) <= EPS_STOP_ratio % if no enough progress, end
        %                 break;
        %             end
        %         end
        % update phi
        [Phi_minus_GP_next,Phi_plus_GP_next] ...
            = Iteration_GP(Adj,Task,InputRate,a_m,b_overhead,StepSize_GP(iter_t),Phi_minus_GP_current,Phi_plus_GP_current,...
            LinkFlow_GP_current,CompFlow_GP_current,t_minus_GP_current,t_plus_GP_current,Delay_deri_GP_current,...
            CompCost_deri_GP_current,Capacity_Saturate_Factor,eps);
        % overwrite
        Phi_minus_GP_current = Phi_minus_GP_next;
        Phi_plus_GP_current = Phi_plus_GP_next;
    end
    
    if Is_Use_SGP % if using scaled gradient projection algorithm
        % update and save current statues
        [LinkFlow_SGP_current,CompFlow_SGP_current,t_minus_SGP_current,t_plus_SGP_current] ...
            = Update_Flow(Adj,M,Task,a_m,b_overhead,InputRate,Phi_minus_SGP_current,Phi_plus_SGP_current,eps);
        
        [Delay_SGP_current,Delay_deri_SGP_current,CompCost_SGP_current,CompCost_deri_SGP_current] ...
            = Update_Cost(Delay_type,CompCost_type,Delay_para,CompCost_para,LinkFlow_SGP_current,CompFlow_SGP_current,eps);
        
        TotalCost_SGP(iter_t) = sum(sum(Delay_SGP_current)) + sum(sum(CompCost_SGP_current));
        %         if iter_t >= 2
        %             if  abs(TotalCost_SGP(iter_t) - TotalCost_SGP(iter_t -1))/TotalCost_SGP(iter_t -1) <= EPS_STOP_ratio % if no enough progress, end
        %                 break;
        %             end
        %         end
        % update phi
        [Phi_minus_SGP_next,Phi_plus_SGP_next] ...
            = Iteration_SGP(Adj,Task,InputRate,a_m,b_overhead,[],Phi_minus_SGP_current,Phi_plus_SGP_current,...
            LinkFlow_SGP_current,CompFlow_SGP_current,t_minus_SGP_current,t_plus_SGP_current,Delay_deri_SGP_current,...
            CompCost_deri_SGP_current,Capacity_Saturate_Factor,eps);
        % overwrite
        Phi_minus_SGP_current = Phi_minus_SGP_next;
        Phi_plus_SGP_current = Phi_plus_SGP_next;
    end
    
    if Is_Use_LB % if using LB, must use GP/SGP to generate initial points and ti-
        if Is_Use_GP
            a_equivalent = zeros(N_node*N_task,1);
    for t_index = 1:N_task
        t_m = Task(t_index,2);
        for node_i = 1:N_node
            a_equivalent((node_i-1)*N_task + t_index) = a_m(t_m) + b_overhead / ( t_minus_GP_current((node_i-1)*N_task + t_index) + eps);
        end
    end
            [TotalCost_Opt_Offline, f_minus_Opt_offline, f_plus_Opt_offline,g_Opt_offline] = ...
        Offline_Optimization(Adj,Task,InputRate,a_equivalent,Delay_type,CompCost_type,Delay_para,CompCost_para,...
        LinkCap,CompCap,f_minus_init,f_plus_init,g_comp_init,Capacity_Saturate_Factor,eps);
    TotalCost_LB(iter_t) = TotalCost_Opt_Offline;
        elseif Is_Use_SGP
            
        else
            error('ERROR: LB must come along with GP/SGP!\n')
        end
    end
    
end

%% plot
figure(2)
if Is_Use_GP
    plot(TotalCost_GP,'b-')
    hold on
end
if Is_Use_SGP
    plot(TotalCost_SGP,'r-')
    hold on
end
if Is_Use_LB
    plot(TotalCost_LB,'y-')
    hold on
end
if Is_Use_Offline
    plot(TotalCost_Offline,'m')
    hold on
end
%hold off

%% Save data and analyze
SaveDataFileName = 'RouteComputeNetwork_Data';
save(SaveDataFileName);

%% test and plot

if Is_debugTest == 1
    
    
    [LinkFlow,CompFlow,t_minus,t_plus] = Update_Flow(Adj,M,Task,a_m,b_overhead,InputRate,Phi_minus_init,Phi_plus_init,eps);
    
    
    
    Task
    InputRate
    %f_minus_init_mat = reshape(f_minus_init,N_node,N_node)'
    %f_plus_init_mat = reshape(f_plus_init,N_node,N_node)'
    %g_comp_init
    if N_task == 1
        Phi_minus_init_mat = reshape(Phi_minus_init,N_node+1,N_node)'
        Phi_plus_init_mat = reshape(Phi_plus_init,N_node,N_node)'
    end
    LinkFlow
    CompFlow
    t_minus
    t_plus
    [Delay,Delay_deri,CompCost,CompCost_deri] = Update_Cost(Delay_type,CompCost_type,Delay_para,CompCost_para,LinkFlow,CompFlow,eps)
    TotalCost_init = sum(sum(Delay)) + sum(sum(CompCost))
    
    a_equivalent = 0.5 * ones(N_node*N_task,1);
    [TotalCost_Opt_Offline, f_minus_Opt_offline, f_plus_Opt_offline,g_Opt_offline] = ...
        Offline_Optimization(Adj,Task,InputRate,a_equivalent,Delay_type,CompCost_type,Delay_para,CompCost_para,...
        LinkCap,CompCap,f_minus_init,f_plus_init,g_comp_init,Capacity_Saturate_Factor,eps);
    
    [Margin_node_minus,Margin_node_plus,Margin_link_minus,Margin_link_plus] =...
        Broadcast_Margin(Adj,Task,InputRate,a_m,b_overhead,Phi_minus_init,Phi_plus_init,Delay_deri,CompCost_deri,t_minus,eps);
    Margin_node_minus
    Margin_node_plus
    if N_task == 1
        Margin_link_minus_mat = reshape(Margin_link_minus,N_node+1,N_node)'
        Margin_link_plus_mat = reshape(Margin_link_plus,N_node,N_node)'
    end
    
    [Is_Blocked_minus,Is_Blocked_plus] = Update_Blocked(Adj,Task,Phi_minus_init,Phi_plus_init,Margin_node_minus,Margin_node_plus,eps);
    if N_task == 1
        Is_Blocked_minus_mat = reshape(Is_Blocked_minus,N_node,N_node)'
        Is_Blocked_plus_mat = reshape(Is_Blocked_plus,N_node,N_node)'
    end
    
    stepsize_test = 0.5;
    [Scale_Mat_minus,Scale_Mat_plus] ... % input 'GP' as type parameter
        = Update_Scale_Matrix(Adj,Task,stepsize_test,Phi_minus_init,Phi_plus_init,Is_Blocked_minus,Is_Blocked_plus, ...
        t_minus,t_plus,Margin_link_minus,Margin_link_plus,'GP',Capacity_Saturate_Factor,eps);
    if N_task == 1
        Scale_Mat_minus_mat_GP = reshape(Scale_Mat_minus,N_node+1,N_node)'
        Scale_Mat_plus_mat_GP = reshape(Scale_Mat_plus,N_node,N_node)'
    end
    %
    % [Scale_Mat_minus,Scale_Mat_plus] ... % input 'GP' as type parameter
    %     = Update_Scale_Matrix(Adj,Task,1,Phi_minus_init,Phi_plus_init,Is_Blocked_minus,Is_Blocked_plus, ...
    %     t_minus,t_plus,Margin_link_minus,Margin_link_plus,'SGP',Capacity_Saturate_Factor);
    % if N_task == 1
    % Scale_Mat_minus_mat_SGP = reshape(Scale_Mat_minus,N_node+1,N_node)'
    % Scale_Mat_plus_mat_SGP = reshape(Scale_Mat_plus,N_node,N_node)'
    % end
    
    [Phi_minus_next,Phi_plus_next] ...
        = Update_Phi_SGP(Adj,Task,Phi_minus_init,Phi_plus_init,Margin_link_minus,Margin_link_plus, ...
        Scale_Mat_minus,Scale_Mat_plus,Is_Blocked_minus,Is_Blocked_plus,eps);
    if N_task == 1
        Phi_minus_next_mat = reshape(Phi_minus_next,N_node+1,N_node)'
        Phi_plus_next_mat = reshape(Phi_plus_next,N_node,N_node)'
    end
    
    [LinkFlow_next,CompFlow_next,t_minus_next,t_plus_next] = Update_Flow(Adj,M,Task,a_m,b_overhead,InputRate,Phi_minus_next,Phi_plus_next,eps)
    [Delay_next,Delay_deri_next,CompCost_next,CompCost_deri_next] = Update_Cost(Delay_type,CompCost_type,Delay_para,CompCost_para,LinkFlow_next,CompFlow_next,eps)
    Cost = sum(sum(Delay)) + sum(sum(CompCost))
    TotalCost_Opt_Offline
    Cost_next = sum(sum(Delay_next)) + sum(sum(CompCost_next))
    
end

%% Appendix: Functions
function [Phi_minus_GP_next,Phi_plus_GP_next] = ...
    Iteration_GP(Adj,Task,InputRate,a_m,b_overhead,StepSize_GP_current,Phi_minus_GP_current,Phi_plus_GP_current,LinkFlow_GP_current,...
    CompFlow_GP_current,t_minus_GP_current,t_plus_GP_current,Delay_deri_GP_current,CompCost_deri_GP_current,Capacity_Saturate_Factor,eps)
% Main iteration for Gradient Projection variable updating

% update marginal delay pD/pr and pD/pt via broadcasting
[Margin_node_minus_GP_current,Mgarin_node_plus_GP_current,Margin_link_minus_GP_current,Margin_link_plus_GP_current] ...
    = Broadcast_Margin(Adj,Task,InputRate,a_m,b_overhead,Phi_minus_GP_current,Phi_plus_GP_current,...
    Delay_deri_GP_current,CompCost_deri_GP_current,t_minus_GP_current,eps);

% calculate bloacked nodes
[Is_Blocked_minus_GP_current,Is_Blocked_plus_GP_current] ...
    = Update_Blocked(Adj,Task,Phi_minus_GP_current,Phi_plus_GP_current,Margin_node_minus_GP_current,Mgarin_node_plus_GP_current,eps);

% calculate scale matrix
[Scale_Mat_minus,Scale_Mat_plus] ... % input 'GP' as type parameter
    = Update_Scale_Matrix(Adj,Task,StepSize_GP_current,Phi_minus_GP_current,Phi_plus_GP_current,Is_Blocked_minus_GP_current,Is_Blocked_plus_GP_current, ...
    t_minus_GP_current,t_plus_GP_current,Margin_link_minus_GP_current,Margin_link_plus_GP_current,'GP',Capacity_Saturate_Factor,eps);

% new phi after projection
[Phi_minus_GP_next,Phi_plus_GP_next] ...
    = Update_Phi_SGP(Adj,Task,Phi_minus_GP_current,Phi_plus_GP_current,Margin_link_minus_GP_current,Margin_link_plus_GP_current, ...
    Scale_Mat_minus,Scale_Mat_plus,Is_Blocked_minus_GP_current,Is_Blocked_plus_GP_current,eps);

end

function [Phi_minus_SGP_next,Phi_plus_SGP_next] = ...
    Iteration_SGP(Adj,Task,InputRate,a_m,b_overhead,~,Phi_minus_SGP_current,Phi_plus_SGP_current,LinkFlow_SGP_current,CompFlow_SGP_current,...
    t_minus_SGP_current,t_plus_SGP_current,Delay_deri_SGP_current,CompCost_deri_SGP_current,Capacity_Saturate_Factor,eps)
% Main iteration for Gradient Projection variable updating

% update marginal delay pD/pr and pD/pt via broadcasting
[Margin_node_minus_SGP_current,Mgarin_node_plus_SGP_current,Margin_link_minus_SGP_current,Margin_link_plus_SGP_current] ...
    = Broadcast_Margin(Adj,Task,InputRate,a_m,b_overhead,Phi_minus_SGP_current,Phi_plus_SGP_current,...
    Delay_deri_SGP_current,CompCost_deri_SGP_current,t_minus_SGP_current,eps);

% calculate bloacked nodes
[Is_Blocked_minus_SGP_current,Is_Blocked_plus_SGP_current] ...
    = Update_Blocked(Adj,Task,Phi_minus_SGP_current,Phi_plus_SGP_current,Margin_node_minus_SGP_current,Mgarin_node_plus_SGP_current,eps);

% calculate scale matrix
[Scale_Mat_minus,Scale_Mat_plus] ... % input 'GP' as type parameter
    = Update_Scale_Matrix(Adj,Task,[],Phi_minus_SGP_current,Phi_plus_SGP_current,Is_Blocked_minus_SGP_current,Is_Blocked_plus_SGP_current, ...
    t_minus_SGP_current,t_plus_SGP_current,Margin_link_minus_SGP_current,Margin_link_plus_SGP_current,'SGP',Capacity_Saturate_Factor,eps);

% new phi after projection
[Phi_minus_SGP_next,Phi_plus_SGP_next] ...
    = Update_Phi_SGP(Adj,Task,Phi_minus_SGP_current,Phi_plus_SGP_current,Margin_link_minus_SGP_current,Margin_link_plus_SGP_current, ...
    Scale_Mat_minus,Scale_Mat_plus,Is_Blocked_minus_SGP_current,Is_Blocked_plus_SGP_current,eps);

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
        for node_j = find(Adj(Adj(node_i,:)))
            Aij_D0_plus(node_i,node_j) = 2/(1-Link_Saturate_Factor);
            if Aij_D0_plus(node_i,node_j) > A_D0
                A_D0 = Aij_D0_plus(node_i,node_j);
            end
        end
    end
    for node_i = 1:N_node
        for node_j = [0 find(Adj(Adj(node_i,:)))]
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
                    Scale_Mat_plus(M_plus_pos) = t_plus(t_plus_pos) / 2 ...
                        * ( Aij_D0_plus(node_i,node_j) +  A_D0 * Unblocked_size_plus * h_plus(node_j));
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
                        Scale_Mat_minus(M_minus_pos) = t_minus(t_minus_pos) / 2 ...
                            * ( Aij_D0_minus(node_i,node_j+1) +  A_D0 * Unblocked_size_minus * h_minus(node_j));
                    end
                end
                
            end
        end
    end
else
    error('ERROR: Wrong Scale Type');
end

end

function [Margin_node_minus,Margin_node_plus,Margin_link_minus,Margin_link_plus] =...
    Broadcast_Margin(Adj,Task,InputRate,a_m,b,Phi_minus,Phi_plus,Delay_deri,CompCost_deri,t_minus,eps)
% Calculate the marginal costs pD/pr and pD/pt+, simulating broadcast
% Phi_minus: column vecto (N*(N+1)*N_task x 1)
% Phi_plus: column vecto (N*N*N_task x 1)
% Delay_deri, CompCost_deri: matrix (N x N)
% Margin_node_minus, Mgarin_node_plus: matrix (N x N_task)
% Margin_link_minus: column vector ( N*(N+1)*N_task x 1 )
% Margin_link_plus: column vector ( N*N*N_task x 1 )
% coeff_t = 0.1;
% offset_t = max(sum(InputRate,2));
coeff_t = 1;
offset_t = 0.01;

N_node = length(Adj);
N_task = size(Task,1);
length_Phi_minus = N_node * (N_node + 1) * N_task;
length_Phi_plus = N_node * N_node * N_task;

% first the plus marginal
Margin_node_plus = zeros(N_node,N_task);
%fprintf('Start: Marginal_plus\n');
for t_index = 1:N_task
    % for each task, using BFS-like starting from destination
    t_dest = Task(t_index,1);
    Is_traversed = zeros(1,N_node); % Marks for completed nodes
    Is_taverse_next = zeros(1,N_node); % Marks for nodes that will determin pD/pt+ next iteration
    DAG_task_plus = zeros(N_node,N_node); % Adj matrix for directed graph of current task
    
    % construncting DAG
    for node_i = 1:N_node
        for node_j = 1:N_node
            Phi_plus_pos = (node_i -1)* N_node * N_task + (node_j-1) * N_task + t_index;
            if Phi_plus(Phi_plus_pos) > eps % mark the link with phi > 0
                if Adj(node_i,node_j) == 0
                    error('ERROR: Positive phi+ on non-link (%d,%d)!\n',node_i,node_j)
                end
                DAG_task_plus(node_i,node_j) = 1;
            end
        end
    end
    
    % starting with destination, traverse nodes with all out-neighbor been traversed
    Is_taverse_next(t_dest) = 1;
    while sum(Is_traversed) < N_node % continue if there is node left
        if sum(Is_taverse_next) == 0 % if not finished but no nodes to traverse next
            error('ERROR: Broadcast abortion!\n')
        end
        %Is_taverse_next
        for node_i = find(Is_taverse_next) % for all nodes to traverse
            for node_j = find(DAG_task_plus(node_i,:)) % update pD/pt+i and delta_ij+
                
                if Is_traversed(node_j) == 0
                    error('ERROR: Traversing node is not ready!\n');
                end
                Phi_plus_pos = (node_i -1)* N_node * N_task + (node_j-1) * N_task + t_index;
                
                % pD/pt+i = \sum_j phiij+ (D'ij + pD/pt+j)
                Margin_node_plus(node_i,t_index) = Margin_node_plus(node_i,t_index) + ...
                    Phi_plus(Phi_plus_pos) * (Margin_node_plus(node_j,t_index) + Delay_deri(node_i,node_j));
                
                % delta_ij+ = D'ij + pD/pt+j
                %Margin_link_plus(Phi_plus_pos) = Margin_node_plus(node_j,t_index) + Delay_deri(node_i,node_j);
            end
            
            Is_traversed(node_i) = 1;
        end
        
        % then mark all nodes to be visit next iter: nodes not traversed but with out-neighbors all traversed.
        Is_taverse_next = zeros(1,N_node);
        for node_i = find(Is_traversed == 0)
            Is_will_traverse = 1;
            for node_j = find(DAG_task_plus(node_i,:)) % for all j: (i,j) in DAG
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

% then the minus marginal
Margin_node_minus = zeros(N_node,N_task);
%fprintf('Start: Marginal_minus\n');
for t_index = 1:N_task
    t_m = Task(t_index,2);
    DAG_task_minus = zeros(N_node,N_node); % Adj matrix for directed graph of current task (minus)
    
    % construncting DAG
    for node_i = 1:N_node
        for node_j = 1:N_node
            Phi_minus_pos = (node_i -1)* (1+N_node) * N_task + node_j * N_task + t_index;
            if Phi_minus(Phi_minus_pos) > eps % mark the link with phi > 0
                if Adj(node_i,node_j) == 0
                    error('ERROR: Positive phi- on non-link (%d,%d)!\n',node_i,node_j)
                end
                DAG_task_minus(node_i,node_j) = 1;
            end
        end
    end
    
    Is_traversed = zeros(1,N_node); % Marks for completed nodes
    Is_taverse_next = zeros(1,N_node);
    % find the start nodes: only forward to local CPU (i.e. no out-nieghbor)
    for node_i = 1:N_node
        if sum(DAG_task_minus(node_i,:)) == 0
            Is_taverse_next(node_i) = 1;
        end
    end
    
    % traverse
    while sum(Is_traversed) < N_node % continue if there is node left
        if sum(Is_taverse_next) == 0 % if not finished but no nodes to traverse next
            error('ERROR: Broadcast abortion!\n')
        end
        %Is_taverse_next
        for node_i = find(Is_taverse_next) % for all nodes to traverse
            % update pD/pr and delta_ij-
            Phi_i0_pos = (node_i-1)*(1+N_node)*N_task + t_index;
            t_minus_pos = (node_i-1)*N_task + t_index;
            
            % phii0*(C'im + (am+b/t-i) pD/pt+i)
            % Note: t-i is multiplied by coeff_t and added by offset_t to prevent infinity value and maintain practical
            Margin_node_minus(node_i,t_index) = Phi_minus(Phi_i0_pos) ...
                * (CompCost_deri(node_i,t_m) + (a_m(t_m)+ b / (coeff_t * t_minus(t_minus_pos) + offset_t)) * Margin_node_plus(node_i,t_index));
            %Margin_link_minus(Phi_i0_pos) ...
            %   = CompCost_deri(node_i,t_m) + (a_m(t_m)+ b / (coeff_t * t_minus(t_minus_pos) + offset_t)) * Margin_node_plus(node_i,t_index);
            
            for node_j = find(DAG_task_minus(node_i,:))
                
                if Is_traversed(node_j) == 0
                    error('ERROR: Traversing node is not ready!\n');
                end
                Phi_minus_pos = (node_i -1)* (1+N_node) * N_task + node_j * N_task + t_index;
                
                % pD/pri = \sum_jneq0 phiij- (D'ij + pD/prj) + phii0*(C'i + (am+b/t-i) pD/pt+i)
                Margin_node_minus(node_i,t_index) = Margin_node_minus(node_i,t_index) + ...
                    Phi_minus(Phi_minus_pos) * (Margin_node_minus(node_j,t_index) + Delay_deri(node_i,node_j));
                
                % delta_ij+ = D'ij + pD/pt+j
                %Margin_link_minus(Phi_minus_pos) = Margin_node_minus(node_j,t_index) + Delay_deri(node_i,node_j);
            end
            
            Is_traversed(node_i) = 1;
        end
        
        % then mark all nodes to be visit next iter: nodes not traversed but with out-neighbors all traversed.
        Is_taverse_next = zeros(1,N_node);
        for node_i = find(Is_traversed == 0)
            Is_will_traverse = 1;
            for node_j = find(DAG_task_minus(node_i,:)) % for all j: (i,j) in DAG
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

% finally globally compute link/compute marginals
Margin_link_minus = zeros(length_Phi_minus,1);
Margin_link_plus = zeros(length_Phi_plus,1);
for t_index = 1:N_task
    t_m = Task(t_index,2);
    for node_i = 1:N_node
        delta_i0_pos = (node_i-1)*(1+N_node)*N_task + t_index;
        t_minus_pos = (node_i-1)*N_task + t_index;
        Margin_link_minus(delta_i0_pos) = ... % delta_i0 = C'im + am * pD/pt+
            CompCost_deri(node_i,t_m) + (a_m(t_m)+ b / (coeff_t * t_minus(t_minus_pos) + offset_t)) * Margin_node_plus(node_i,t_index);
        for node_j = find(Adj(node_i,:)) % for all j: (i,j) is link
            delta_minus_pos = (node_i-1)*(1+N_node)*N_task + node_j*N_task + t_index;
            delta_plus_pos = (node_i-1)*N_node*N_task + (node_j-1)*N_task + t_index;
            Margin_link_minus(delta_minus_pos) = Delay_deri(node_i,node_j) + Margin_node_minus(node_j,t_index);% delta_ij- = D'ij + pD/prj
            Margin_link_plus(delta_plus_pos) = Delay_deri(node_i,node_j) + Margin_node_plus(node_j,t_index);% delta_ij+ = D'ij + pD/pt+j
        end
    end
end
end

function [Phi_minus_new,Phi_plus_new] = ...
    Update_Phi_SGP(Adj,Task,Phi_minus,Phi_plus,Margin_link_minus,Margin_link_plus,Scale_Mat_minus,Scale_Mat_plus,Is_Blocked_minus,Is_Blocked_plus,eps)
% Calculate projected variable using (26),(27) in overleaf, or (36) in Yufang paper
% Note: to avoid singularity or degenerated objective, could add a positive offset to scaling matrix
% Phi_minus_new: column vector (N*(N+1)*N_task x 1);
% Phi_plus_new: column vector (N*N*N_task x 1);
options =  optimoptions(@quadprog,'Display','off');
Offset_ratio = 0.01; % ratio of offset to the maximum entry
Offset_fix = 0.01; % fixed offset

N_node = length(Adj);
N_task = size(Task,1);
Phi_minus_new = zeros(size(Phi_minus));
Phi_plus_new = zeros(size(Phi_plus));

for t_index = 1:N_task
    t_dest = Task(t_index,1);
    for node_i = 1:N_node
        % solve convex quadratic optimization:
        % min  deltaT*(x-phit) + (x-phit)T*M/2*(x-phit), subject to xij >= 0 and sumj(j not blocked) xij = 1, sumj(j blocked)xij = 0
        % method: x = quadprog(H,f,A,b,Aeq,beq,lb,ub): min_x 1/2*xT*H*x + fT*x, such that Ax<=b, Aeq*x=beq, lb<=x<=ub
        % substitute: H = M, f = delta_i, and s.t. xj>=-phitj, sumj(j blocked) xj = - sumj(j blocked) phitj,
        % sumj(j not blocked)xj = - sumj(j not blocked)phitj + 1
        
        % first plus
        delta_plus_pos_list = ((node_i-1)*N_node*N_task+t_index) :N_task: ((node_i-1)*N_node*N_task+(N_node-1)*N_task+t_index);
        delta_plus = Margin_link_plus(delta_plus_pos_list);
        M_plus_pos_list = delta_plus_pos_list;
        M_diag_plus = Scale_Mat_plus(M_plus_pos_list);
        Diag_max = max(M_diag_plus);
        Offset_plus = Diag_max * Offset_ratio + Offset_fix;
        M_plus = diag(M_diag_plus + Offset_plus);
        A_eq_plus = zeros(2,N_node); % firts row is blocked, second row unblocked
        if t_dest == node_i
            b_eq_plus = [0;0]; % if is the destination, no out flow
        else
            b_eq_plus = [0;1];
        end
        lb_plus = zeros(N_node,1);
        for node_j = 1:N_node
            Phi_plus_pos = (node_i-1)*N_node*N_task + (node_j-1)*N_task + t_index;
            Is_blocked_plus_pos = Phi_plus_pos;
            lb_plus(node_j) = -1 * Phi_plus(Phi_plus_pos);
            if Is_Blocked_plus(Is_blocked_plus_pos) == 1 % if no link ij or j is blocked
                A_eq_plus(1,node_j) = 1;
                b_eq_plus(1) = b_eq_plus(1) - Phi_plus(Phi_plus_pos);
            else % if ij is feasible out-link
                A_eq_plus(2,node_j) = 1;
                b_eq_plus(2) = b_eq_plus(2) - Phi_plus(Phi_plus_pos);
            end
        end
        [phi_plus_shift_idm,fval,exitflag_1] = quadprog(M_plus,delta_plus,[],[],A_eq_plus,b_eq_plus,lb_plus,[],[],options);
        if (exitflag_1 <0)
            error('ERROR: Quadratic progamming fail!\n')
        end
        phi_plus_new_idm = phi_plus_shift_idm - lb_plus; % phi_new = x + phi_old
        
        % then minus. Note: inlcude i=0
        delta_minus_pos_list = ((node_i-1)*(1+N_node)*N_task+t_index) :N_task: ((node_i-1)*(1+N_node)*N_task+N_node*N_task+t_index);
        delta_minus = Margin_link_minus(delta_minus_pos_list);
        M_minus_pos_list = delta_minus_pos_list;
        M_diag_minus = Scale_Mat_minus(M_minus_pos_list);
        Diag_max = max(M_diag_minus);
        Offset_minus = Diag_max * Offset_ratio + Offset_fix;
        M_minus = diag(M_diag_minus + Offset_minus);
        A_eq_minus = zeros(2,1+N_node); % firts row is blocked, second row unblocked
        b_eq_minus = [0;1];
        lb_minus = zeros(1+N_node,1);
        for node_j = 0:N_node % note: start from 0
            Phi_minus_pos = (node_i-1)*(1+N_node)*N_task + node_j*N_task + t_index;
            lb_minus(node_j+1) = -1 * Phi_minus(Phi_minus_pos);
            
            if node_j == 0 % for j=0, not blocked
                A_eq_minus(2,node_j+1) = 1;
                b_eq_minus(2) = b_eq_minus(2) - Phi_minus(Phi_minus_pos);
            else
                Is_blocked_minus_pos = (node_i-1)*N_node*N_task + (node_j-1)*N_task + t_index;
                if Is_Blocked_minus(Is_blocked_minus_pos) == 1
                    A_eq_minus(1,node_j+1) = 1;
                    b_eq_minus(1) = b_eq_minus(1) - Phi_minus(Phi_minus_pos);
                else
                    A_eq_minus(2,node_j+1) = 1;
                    b_eq_minus(2) = b_eq_minus(2) - Phi_minus(Phi_minus_pos);
                end
            end
        end
        [phi_minus_shift_idm,fval,exitflag_2] = quadprog(M_minus,delta_minus,[],[],A_eq_minus,b_eq_minus,lb_minus,[],[],options);
        if (exitflag_2 <0)
            error('ERROR: Quadratic progamming fail!\n')
        end
        phi_minus_new_idm = phi_minus_shift_idm - lb_minus;
        
        
        % assign entries
        for node_j = 0:N_node
            if node_j == 0
                Phi_i0_pos = (node_i -1)*(1+N_node)*N_task + t_index;
                Phi_minus_new(Phi_i0_pos) = phi_minus_new_idm(1);
            else
                Phi_minus_pos = (node_i -1)*(1+N_node)*N_task + node_j*N_task + t_index;
                Phi_minus_new(Phi_minus_pos) = phi_minus_new_idm(node_j +1 );
                Phi_plus_pos = (node_i -1)*N_node*N_task + (node_j-1)*N_task + t_index;
                Phi_plus_new(Phi_plus_pos) = phi_plus_new_idm(node_j);
            end
        end
    end
end

end

function [Is_Blocked_minus,Is_Blocked_plus] = ...
    Update_Blocked(Adj,Task,Phi_minus,Phi_plus,Margin_node_minus,Margin_node_plus,eps)
% Calculate the blocked nodes accroding to pD/pr and pD/pt+.
% Blocked_Mat_minus, Blocked_Mat_plus: column vector (N*N*N_task x 1), (i,j,t) denotes if j is blocked to i wrt task t.
% Method:   First calculate improper nodes by BFS (i.e., containing an improper path to destination/sink)
%           Then assign entries to blocked_mat: either no ij link, or phiij = 0 and j is improper.
N_node = length(Adj);
N_task = size(Task,1);

length_Is_Blocked_minus = N_node * N_node * N_task;
length_Is_Blocked_plus = N_node * N_node * N_task;
Is_Blocked_minus = zeros(length_Is_Blocked_minus,1);
Is_Blocked_plus = zeros(length_Is_Blocked_plus,1);
for t_index = 1:N_task % separate for each task
    Is_improper_plus = zeros(1,N_node); % improper nodes indicator
    Is_improper_minus = zeros(1,N_node);
    Is_improper_link_plus = zeros(N_node,N_node); % improper links indicator
    Is_improper_link_minus = zeros(N_node,N_node);
    
    % first mark all improper links
    for node_i = 1:N_node
        for node_j = find(Adj(node_i,:)) % all links i,j
            Phi_plus_pos = (node_i-1)*N_node*N_task + (node_j-1)*N_task + t_index;
            Phi_minus_pos = (node_i-1)*(1+N_node)*N_task + node_j*N_task + t_index;
            if (Phi_plus(Phi_plus_pos) > eps) ... % first plus: if phiij+ >0 and pD/pt+i < pD/pt+j
                    && (Margin_node_plus(node_j,t_index) - Margin_node_plus(node_i,t_index) > eps)
                Is_improper_link_plus(node_i,node_j) = 1;
            end
            if (Phi_minus(Phi_minus_pos) > eps) ... % then minus: if phiij- >0 and pD/pri < pD/prj
                    && (Margin_node_minus(node_j,t_index) - Margin_node_minus(node_i,t_index) > eps)
                Is_improper_link_minus(node_i,node_j) = 1;
            end
        end
    end
    
    % then BSF (flood for at most N iterations)
    for iter = 1:N_node
        Is_changed = 0;
        % plus
        for node_i = find(Is_improper_plus == 0) % exam all nodes not marked as improper,
            for node_j = find(Adj(node_i,:))
                if Is_improper_link_plus(node_i,node_j) == 1 % mark if it has improper out-link
                    Is_improper_plus(node_i) = 1;
                    Is_changed = 1;
                    break;
                elseif Is_improper_plus(node_j) == 1% mark if it has improper out-neighbor
                    Is_improper_plus(node_i) = 1;
                    Is_changed = 1;
                    break;
                end
            end
        end
        
        % minus
        for node_i = find(Is_improper_minus == 0) % exam all nodes not marked as improper,
            for node_j = find(Adj(node_i,:)) % all j: (i,j) is link
                if Is_improper_link_minus(node_i,node_j) == 1 % mark if it has improper out-link
                    Is_improper_minus(node_i) = 1;
                    Is_changed = 1;
                    break;
                elseif Is_improper_minus(node_j) == 1% mark if it has improper out-neighbor
                    Is_improper_minus(node_i) = 1;
                    Is_changed = 1;
                    break;
                end
            end
        end
        if Is_changed == 0
            break;
        end
    end
    
    % then assign to output
    for node_i = 1:N_node
        for node_j = 1:N_node
            Is_Blocked_pos = (node_i-1)*N_node*N_task + (node_j-1)*N_task + t_index;
            if Adj(node_i,node_j) == 0 % if not a link, block both - and +
                Is_Blocked_minus(Is_Blocked_pos) = 1;
                Is_Blocked_plus(Is_Blocked_pos) = 1;
            else % if ij is a link
                Phi_minus_pos = (node_i -1)*(1+N_node)*N_task + node_j*N_task + t_index;
                Phi_plus_pos = (node_i-1)*N_node*N_task + (node_j -1)*N_task + t_index;
                if Phi_minus(Phi_minus_pos) < eps % only block nodes with no current flow
                    if Margin_node_minus(node_j,t_index) - Margin_node_minus(node_i,t_index) > eps % if there is link ij, but j has higher marginal
                        Is_Blocked_minus(Is_Blocked_pos) = 1;
                    end
                    if Is_improper_minus(node_j) == 1
                        Is_Blocked_minus(Is_Blocked_pos) = 1;
                    end
                end
                
                if Phi_plus(Phi_plus_pos) < eps
                    if Margin_node_plus(node_j,t_index) - Margin_node_plus(node_i,t_index) > eps % if there is link ij, but j has higher marginal
                        Is_Blocked_plus(Is_Blocked_pos) = 1;
                    end
                    if Is_improper_plus(node_j) == 1 % if there is link ij, but j is improper, block ij
                        Is_Blocked_plus(Is_Blocked_pos) = 1;
                    end
                end
                
                
            end
        end
    end
end
end

function [Task,InputRate] = ...
    Task_Generator_random(Adj,N_task,M,rate_min,rate_max,inpute_nodes_number_max)
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

function [Phi_minus,Phi_plus] = ...
    FtoPhi(Adj,Task,InputRate,a_m,b,f_minus,f_plus,g_comp,eps)

N_node = length(Adj);
N_task = size(Task,1);
% first calculate t- ,t+ using f-,f+,g
Mat_FtoT = kron(ones(1,N_node),eye(N_node*N_task));
t_minus = Mat_FtoT * f_minus + reshape(InputRate',[],1); % ti- = sum_j fji- + ri
extra_flow = zeros(size(g_comp)); % \alpha(gidm)
for node_i = 1:N_node
    for t_index = 1:N_task
        g_pos = (node_i-1)*N_task + t_index;
        if g_comp(g_pos) <= eps
            extra_flow(g_pos) = 0;
        else
            extra_flow(g_pos) = g_comp(g_pos) * a_m(Task(t_index,2)) + b;
        end
    end
end
t_plus =  Mat_FtoT * f_plus + extra_flow;

% then calculate phi
length_Phi_minus = N_node * (1+N_node) * N_task; % column vector odered as i ,j ,(d,m) and j can be 0
length_Phi_plus = N_node * N_node * N_task; % column vector
Phi_minus = zeros(length_Phi_minus,1);
Phi_plus = zeros(length_Phi_plus,1);
for t_index = 1:N_task
    t_dest = Task(t_index,1);
    for node_i = 1:N_node
        t_pos = (node_i-1)*N_task + t_index;
        
        % first assign phi-
        if t_minus(t_pos) < eps % if t-idm = 0, assign all phi-ij to j=0
            Phi_i0_pos = (node_i-1)* (N_node+1) * N_task + t_index;
            Phi_minus(Phi_i0_pos) = 1;
        else % if t-idm > 0, assign regularly
            g_comp_pos = (node_i-1)*N_task + t_index;
            Phi_i0_pos = (node_i-1)* (N_node+1) * N_task + t_index;
            Phi_minus(Phi_i0_pos) = g_comp(g_comp_pos) / t_minus(t_pos);
            for node_j = 1:N_node
                f_pos = (node_i -1)*N_node*N_task + (node_j-1)*N_task + t_index;
                Phi_minus_pos = (node_i -1)*(N_node+1)*N_task + node_j*N_task + t_index;
                Phi_minus(Phi_minus_pos) = f_minus(f_pos) / t_minus(t_pos);
            end
        end
        
        % then assign phi+
        if node_i ~= t_dest % if i is the destination, do nothing (keep all-0)
            if t_plus(t_pos) < eps % if no flow, pick some out-node to route
                out_node_list = find(Adj(node_i,:)); % list of out-nodes
                [val,out_node_pos] = min(abs(out_node_list - t_dest));% pick node j that has closed index number with destination, to prevent loop
                out_node = out_node_list(out_node_pos);
                Phi_plus_pos = (node_i -1)*N_node*N_task + (out_node-1)*N_task + t_index;
                Phi_plus(Phi_plus_pos) = 1;
            else % if i is not the destination
                for node_j = 1:N_node
                    f_pos = (node_i -1)*N_node*N_task + (node_j-1)*N_task + t_index;
                    Phi_plus_pos = (node_i -1)*N_node*N_task + (node_j-1)*N_task + t_index;
                    Phi_plus(Phi_plus_pos) = f_plus(f_pos) / t_plus(t_pos);
                end
            end
        end
    end
end
end

function [Is_Success,f_minus_init,f_plus_init,g_comp_init] = ...
    Init_Generator_MILP(Adj,Task,InputRate,a_m,b,LinkCap,CompCap,eps)
% generate a vaild initial state by MILP on flow domain, with random
% objective
N_node = length(Adj);
N_task = size(Task,1);
length_f_minus = N_node * N_node * N_task; % column vector ordered as i, j, (d,m)
length_f_plus = N_node * N_node * N_task; % column vector
length_g = N_node * N_task; % column vector ordered as i, (d,m)

UB = sum(sum(InputRate))*1.5; % upper bound for all possible computation flow
computation_cap = UB; %(optional) a sharp upperbound additional to the linear computation cost

length_x_MILP = N_node^2 * 2 * N_task + N_node * N_task * 2; % dimension of x vetcor, arranged with f-ijt,f+ijt,git,wit

% min fx, st. Ax<=b

f_MILP = ones(1,length_x_MILP); % random objective
%f_MILP = [1 zeros(1,length_x_MILP-1)];
%f_MILP = rand(1,length_x_MILP);

length_constraints = N_node * N_task * 3 + N_node * N_node + N_node; % flow conserv for f- at i,f+ at i, g<U*w at i. And LinkCap for all (i,j), CompCap for all i.
A_MILP = zeros(length_constraints,length_x_MILP);
b_MILP = zeros(length_constraints,1);
for node_i = 1:N_node
    for t = 1:N_task
        index_conserv_minus = (node_i-1)*N_task + t;
        index_conserv_plus = N_node * N_task + (node_i-1)*N_task + t;
        index_gleqw = 2*N_node*N_task + (node_i-1)*N_task + t;
        
        % first set the conserv for f-
        b_MILP(index_conserv_minus) = -1* InputRate(node_i,t);
        index_g = N_node^2 * 2 * N_task + (node_i-1)*N_task + t;
        A_MILP(index_conserv_minus,index_g) = -1;
        for node_j = 1:N_node
            if Adj(node_j,node_i) > 0 % in-going links (j,i)
                index_f_minus = (node_j-1)*N_node*N_task + (node_i-1)*N_task + t;
                A_MILP(index_conserv_minus,index_f_minus) = 1;
            end
            if Adj(node_i,node_j) > 0 % out-going links (i,j)
                index_f_minus = (node_i-1)*N_node*N_task + (node_j-1)*N_task + t;
                A_MILP(index_conserv_minus,index_f_minus) = -1;
            end
        end
        
        % then the conserv for f+
        if Task(t,1) ~= node_i % \sum_j fij+ >= \sum_j fji+ + am*gi + b*wi
            index_g = N_node^2 * 2 * N_task + (node_i-1)*N_task + t;
            A_MILP(index_conserv_plus,index_g) = a_m(Task(t,2));
            index_w = N_node^2 * 2 * N_task + N_node * N_task + (node_i-1)*N_task + t;
            A_MILP(index_conserv_plus,index_w) = b;
            for node_j = 1:N_node
                if Adj(node_j,node_i) > 0 % in-going links
                    index_f_plus = N_node^2*N_task + (node_j-1)*N_node*N_task + (node_i-1)*N_task + t;
                    A_MILP(index_conserv_plus,index_f_plus) = 1;
                end
            end
        end
        for node_j = 1:N_node
            if Adj(node_i,node_j) >0 % out-going links
                index_f_plus = N_node^2*N_task + (node_i-1)*N_node*N_task + (node_j-1)*N_task + t;
                A_MILP(index_conserv_plus,index_f_plus) = -1;
            end
        end
        
        % then the constraint gi <= UB * wi
        index_g = N_node^2 * 2 * N_task + (node_i-1)*N_task + t;
        index_w = N_node^2 * 2 * N_task + N_node * N_task + (node_i-1)*N_task + t;
        A_MILP(index_gleqw,index_g) = 1;
        A_MILP(index_gleqw,index_w) = -1 * UB;
    end
    
    % then the link capacity constraint for all (i,j), sum_t (f-ijt + f+ijt) <= LinkCap
    for node_j = 1:N_node
        index_LinkCap = N_node * N_task * 3 + (node_i -1)*N_node + node_j;
        for t = 1:N_task
            index_f_minus =  (node_i - 1)*N_node*N_task + (node_j -1)*N_task + t;
            index_f_plus =  N_node^2*N_task + (node_i - 1)*N_node*N_task + (node_j -1)*N_task + t;
            A_MILP(index_LinkCap,index_f_minus) = 1;
            A_MILP(index_LinkCap,index_f_plus) = 1;
        end
        b_MILP(index_LinkCap) = LinkCap(node_i,node_j);
    end
    
    % then the computation capacity for all nodes, sum_t git <= CompCap
    index_CompCap = N_node * N_task *3 + N_node * N_node + node_i;
    for t = 1:N_task
        index_g_it = N_node^2*N_task*2 + (node_i -1)*N_task + t;
        A_MILP(index_CompCap,index_g_it) = 1;
    end
    b_MILP(index_CompCap) = CompCap(node_i);
end

x_lb_MILP = zeros(1,length_x_MILP);
x_ub_MILP = [UB*ones(1,N_node^2 * 2 * N_task) computation_cap*ones(1,N_node * N_task) ones(1,N_node * N_task)];
int_index = N_node^2 * 2 * N_task + N_node * N_task + 1 : N_node^2 * 2 * N_task + 2 * N_node * N_task;
options =  optimoptions(@intlinprog,'Display','off');
[x_opt_MILP, cost_opt_MILP,exitflag] = intlinprog(f_MILP,int_index,A_MILP,b_MILP,[],[],x_lb_MILP,x_ub_MILP,options);

if exitflag <= 0 % No feasible solution or unbounded
    Is_Success = 0;
    f_minus_init = zeros(length_f_minus,1);
    f_plus_init = zeros(length_f_plus,1);
    g_comp_init = zeros(length_g,1);
else
    Is_Success = 1;
    x_opt_MILP = reshape(x_opt_MILP,[],1);
    f_minus_init = x_opt_MILP(1:length_f_minus);
    f_plus_init = x_opt_MILP(length_f_minus +1 : length_f_minus + length_f_plus);
    g_comp_init = x_opt_MILP(length_f_minus + length_f_plus + 1: length_f_minus + length_f_plus + length_g);
    w_comp_init = x_opt_MILP(length_f_minus + length_f_plus + length_g + 1: end);
end

end

function [LinkFlow,CompFlow,t_minus,t_plus] = ...
    Update_Flow(Adj,M,Task,a_m,b,InputRate,Phi_minus,Phi_plus,eps)
% calculate flow Fij, Gim and Ti+- from routing variable phi
% LinkFlow: matrix (N x N)
% CompFlow: matrix (N x M)
% t_minus, t_plus: column vector (N*N_task, 1)
% Adj: adj matrix; M: computation types; Task: tasks [d,m]; Phi is given in column vector
% Result is computed in order: t- -> g -> t+ -> f- f+ -> F
N_node = length(Adj);
N_task = size(Task,1);
MAX_PATH_LEN = N_node;

length_Phi_minus = N_node * (1+N_node) * N_task; % column vector odered as i ,j ,(d,m) and j can be 0
length_Phi_plus = N_node * N_node * N_task; % column vector
length_f_minus = N_node * N_node * N_task; % column vector ordered as i, j, (d,m)
length_f_plus = N_node * N_node * N_task; % column vector
length_g = N_node * N_task; % column vector ordered as i, (d,m)
length_t_minus = N_node * N_task; % column vector
length_t_plus = N_node * N_task; % column vector

% Method: iterate calculate tidm = \sum j tjdm * phi_ji, start from data source, run for at most N times
t_minus = zeros(length_t_minus,1);
t_plus = zeros(length_t_plus,1);
f_minus = zeros(length_f_minus,1);
f_plus = zeros(length_f_plus,1);
g_comp = zeros(length_g,1);
LinkFlow = zeros(N_node,N_node); % Fij, matrix
CompFlow = zeros(N_node,M); % Gim, matrix

for t_index = 1:N_task % for each task
    t_dest = Task(t_index,1);
    t_m = Task(t_index,2);
    
    % extract phi into matrix
    Phi_minus_t = zeros(N_node,N_node); % phi- for current task, matrix, exclude phi_i0
    Phi_plus_t = zeros(N_node,N_node);
    Phi_0_t = zeros(N_node,1); % phi_i0 for current task, column vector
    for node_i = 1:N_node
        phi_i0_pos = (node_i-1) * (1+N_node)*N_task + t_index;
        Phi_0_t(node_i) = Phi_minus(phi_i0_pos);
        for node_j = 1:N_node
            phi_minus_pos = (node_i-1) * (1+N_node) * N_task + node_j * N_task + t_index; % Note to avoid phi_i0
            phi_plus_pos = (node_i-1) * N_node * N_task + (node_j-1) * N_task + t_index;
            Phi_minus_t(node_i,node_j) = Phi_minus(phi_minus_pos);
            Phi_plus_t(node_i,node_j) = Phi_plus(phi_plus_pos);
        end
    end
    
    % calculate t-
    t_minus_t = InputRate(:,t_index); % t- for current task, column vector, vary every iteration, init with ridm
    for iter = 1:MAX_PATH_LEN
        t_minus_t_new = (t_minus_t' * Phi_minus_t)' + InputRate(:,t_index); % calculate new t using old t: t-j = sum_i t-i * phi_ij + rj
        if sum(abs(t_minus_t_new - t_minus_t)) <= eps % if not changing in one iteration
            break
        end
        t_minus_t = t_minus_t_new;
    end
    
    % calculate g
    g_comp_t = t_minus_t .* Phi_0_t; % gidm, column vector
    
    % calculate t+
    t_plus_t = (a_m(t_m)* g_comp_t + b) .* (g_comp_t > eps); % init with extra flows am*g + b
    for iter = 1:MAX_PATH_LEN
        t_plus_t_new = (t_plus_t' * Phi_plus_t)' + (a_m(t_m)* g_comp_t + b) .* (g_comp_t > eps);
        if sum(abs(t_plus_t_new - t_plus_t)) <= eps
            break
        end
        t_plus_t = t_plus_t_new;
    end
    
    % assign entries
    for node_i = 1:N_node
        g_pos = (node_i-1)*N_task + t_index;
        g_comp(g_pos) = g_comp_t(node_i);
        t_pos = (node_i-1)*N_task + t_index;
        t_minus(t_pos) = t_minus_t(node_i);
        t_plus(t_pos) = t_plus_t(node_i);
        CompFlow(node_i,t_m) = CompFlow(node_i,t_m) + g_comp_t(node_i);
        for node_j = 1:N_node
            f_pos = (node_i-1)*N_node*N_task + (node_j-1)*N_task + t_index;
            f_minus(f_pos) = t_minus_t(node_i) * Phi_minus_t(node_i,node_j);
            f_plus(f_pos) = t_plus_t(node_i) * Phi_plus_t(node_i,node_j);
            LinkFlow(node_i,node_j) = LinkFlow(node_i,node_j) + f_minus(f_pos) + f_plus(f_pos);
        end
    end
end

end

function [Delay,Delay_deri,CompCost,CompCost_deri] = ...
    Update_Cost(Delay_type,CompCost_type,Delay_para,CompCost_para,LinkFlow,CompFlow,eps)
% delay calculation. Link_Flow: N_node * N_node; Delay: N * N; Delay_deri: N * N
% Computation cost calculation, Comp_Flow: N_node * M; CompCost: N * 1; CompCost_deri: N * M
% Queueing Delay: D = F/(C-F)
N_node = length(LinkFlow);
Delay = zeros(size(LinkFlow));
Delay_deri = zeros(size(LinkFlow));
if strcmp(Delay_type,'queue')
    Cap =  Delay_para;
    for node_i = 1:N_node
        for node_j = 1:N_node
            c = Cap(node_i,node_j);
            f = LinkFlow(node_i,node_j);
            if c < eps % if there is no link ij
                if f < eps % if no flow
                    Delay(node_i,node_j) = 0;
                    Delay_deri(node_i,node_j) = 0;
                else
                    Delay(node_i,node_j) = Inf;
                    Delay_deri(node_i,node_j) = Inf;
                    error('ERROR: Exceed Link Capacity');
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
    Delay = Delay_para .* LinkFlow;
    Delay_deri = Delay_para;
end

if strcmp(CompCost_type,'sum_queue')
    Cap =  CompCost_para;
    Comp_Flow_Sum = sum(CompFlow,2);
    if min(min(Cap - Comp_Flow_Sum)) <= 0
        CompCost = Inf * Comp_Flow_Sum;
        CompCost_deri = Inf * CompFlow;
        error('ERROR: Exceed Computation Capacity');
    else
        CompCost = Comp_Flow_Sum ./ (Cap - Comp_Flow_Sum);
        CompCost_deri = Cap ./ (Cap - Comp_Flow_Sum).^2;
        CompCost_deri = kron(CompCost_deri,ones(1,size(CompFlow,2)));
    end
elseif strcmp(CompCost_type,'sum_linear')
    CompCost = CompCost_para .* sum(CompFlow,2);
    CompCost_deri = CompCost_para;
    CompCost_deri = kron(CompCost_deri,ones(1,size(CompFlow,2)));
end

end

function Is_Valid = ...
    Check_Phi_Valid(Adj,Task,InputRate,Phi_minus,Phi_plus)
% check if given phi is valid:
% 1: are flow conervation satisfied
% 2: are + / - flow loop-free

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




