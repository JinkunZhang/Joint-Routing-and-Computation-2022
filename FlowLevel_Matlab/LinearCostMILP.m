% Flow model solver for the joint routing/compuation with b \neq 0 and
% linear cost.
% Test if the solutions of MILP and relaxed LP are equal.
% Created by Jinkun, Feb 08 2022.

%% Setting
N_node = 10; % number of nodes
N_task = 5;
M = 3; % number of computation type
p_extralink = 0.2; % graph is a linear network with extra random link
cost_max = 10; % max cost for link and computation (uniform)
inpute_nodes_number_max = 1; % the maximum input node number of each task
rate_max = 5; % max input rate of task (uniform)
computation_cap = 3; % (optional) a sharp upperbound additional to the linear computation cost

% data-result convert rate
a_m = rand(1,M);
b_overhead = 1; % b in paper

% random adjacent matrix and cost (dij and cim)
Adj_cost = (rand(N_node) < p_extralink);
for i = 1:N_node-1
    Adj_cost(i,i+1) = 1;
end
Adj_cost = Adj_cost .* rand(N_node) * cost_max;
for i = 1:N_node
    Adj_cost(i,i) = 0;
   for j = 1:i-1
       Adj_cost(i,j) = Adj_cost(j,i);
   end
end
Computation_cost = rand(N_node,M) * cost_max; % assume each node and each computation has a different cost

% random tasks
Task = zeros(N_task,2); 
for i = 1:N_task
   Task(i,1) = randi(N_node); % each row of Task is a (d,m)
   Task(i,2) = randi(M);
end

% random input rate
input_nodes_number = randi(inpute_nodes_number_max,1,N_task); % input node number of each task
input_rate = zeros(N_node,N_task); 
for i_task = 1:N_task
   input_n_vec = randperm(N_node, input_nodes_number(i_task));
   for i_node = 1:length(input_n_vec)
       input_rate(input_n_vec(i_node),i_task) = rate_max; %rand() * rate_max;
   end
end

%% generate and plot graph
G = graph(Adj_cost);
plot(G)

%% Solve the flow model MILP
UB = 1e5; % upper bound for all possible computation flow

length_x_MILP = N_node^2 * 2 * N_task + N_node * N_task * 2; % dimension of x vetcor, arranged with f-ijt,f+ijt,git,wit

% min fx, st. Ax<=b

f_MILP = zeros(1,length_x_MILP); % sum_dm [sum_ij dij*(f-ijdm + f+ijdm) + cim*gi] 
for t = 1:N_task % for all (d,m)
    for i = 1:N_node % for all i
      for j = 1:N_node % for all (i,j)
          index_f_minus = (i-1)*N_node*N_task + (j-1)*N_task + t;
          index_f_plus = N_node^2 * N_task + (i-1)*N_node*N_task + (j-1)*N_task + t;
          f_MILP(index_f_minus) = Adj_cost(i,j);
          f_MILP(index_f_plus) = Adj_cost(i,j);
      end
      
      index_g = N_node^2 * 2 * N_task + (i-1)*N_task + t;
      f_MILP(index_g) = Computation_cost(Task(t,2));
    end
end

length_constraints = N_node * N_task * 3; % flow conserv for f- at i,f+ at i, g<U*w at i.
A_MILP = zeros(length_constraints,length_x_MILP);
b_MILP = zeros(length_constraints,1);
for i = 1:N_node
   for t = 1:N_task
       index_conserv_minus = (i-1)*N_task + t;
       index_conserv_plus = N_node * N_task + (i-1)*N_task + t;
       index_gleqw = 2*N_node*N_task + (i-1)*N_task + t;
       
% first set the conserv for f-
       b_MILP(index_conserv_minus) = -1* input_rate(i,t);
       index_g = N_node^2 * 2 * N_task + (i-1)*N_task + t;
       A_MILP(index_conserv_minus,index_g) = -1;
       for j = 1:N_node
           if Adj_cost(j,i) > 0 % in-going links (j,i)
               index_f_minus = (j-1)*N_node*N_task + (i-1)*N_task + t;
               A_MILP(index_conserv_minus,index_f_minus) = 1;
           end
           if Adj_cost(i,j) > 0 % out-going links (i,j)
              index_f_minus = (i-1)*N_node*N_task + (j-1)*N_task + t;
              A_MILP(index_conserv_minus,index_f_minus) = -1;
           end
       end

% then the conserv for f+
        if Task(t,1) ~= i
            index_g = N_node^2 * 2 * N_task + (i-1)*N_task + t;
            A_MILP(index_conserv_plus,index_g) = a_m(Task(t,2));
            index_w = N_node^2 * 2 * N_task + N_node * N_task + (i-1)*N_task + t;
            A_MILP(index_conserv_plus,index_w) = b_overhead;
            for j = 1:N_node
                if Adj_cost(j,i) > 0 % in-going links
                   index_f_plus = N_node^2*N_task + (j-1)*N_node*N_task + (i-1)*N_task + t;
                   A_MILP(index_conserv_plus,index_f_plus) = 1;
                end
            end
        end
        for j = 1:N_node
           if Adj_cost(i,j) >0 % out-going links
               index_f_plus = N_node^2*N_task + (i-1)*N_node*N_task + (j-1)*N_task + t;
               A_MILP(index_conserv_plus,index_f_plus) = -1;
           end
        end

% then the constraint gi <= UB * wi
        index_g = N_node^2 * 2 * N_task + (i-1)*N_task + t;
        index_w = N_node^2 * 2 * N_task + N_node * N_task + (i-1)*N_task + t;
        A_MILP(index_gleqw,index_g) = 1;
        A_MILP(index_gleqw,index_w) = -1 * UB;
   end
end

x_lb_MILP = zeros(1,length_x_MILP);
x_ub_MILP = [UB*ones(1,N_node^2 * 2 * N_task) computation_cap*ones(1,N_node * N_task) ones(1,N_node * N_task)];
int_index = N_node^2 * 2 * N_task + N_node * N_task + 1 : N_node^2 * 2 * N_task + 2 * N_node * N_task;
[x_opt_MILP, cost_opt_MILP] = intlinprog(f_MILP,int_index,A_MILP,b_MILP,[],[],x_lb_MILP,x_ub_MILP);

%% solve the equivalent LP for case when active input node is 1
UB = 1e5;
length_x_LP = N_node^2 * 2 * N_task + N_node * N_task; % dimension of x vetcor, arranged with f-ijt,f+ijt,git

% min fx, st. Ax<=b

f_LP = zeros(1,length_x_LP); % sum_dm [sum_ij dij*(f-ijdm + f+ijdm) + cim*gi] 
for t = 1:N_task % for all (d,m)
    for i = 1:N_node % for all i
      for j = 1:N_node % for all (i,j)
          index_f_minus = (i-1)*N_node*N_task + (j-1)*N_task + t;
          index_f_plus = N_node^2 * N_task + (i-1)*N_node*N_task + (j-1)*N_task + t;
          f_LP(index_f_minus) = Adj_cost(i,j);
          f_LP(index_f_plus) = Adj_cost(i,j);
      end
      
      index_g = N_node^2 * 2 * N_task + (i-1)*N_task + t;

      f_LP(index_g) = Computation_cost(Task(t,2));
    end
end

length_constraints = N_node * N_task * 2; % flow conserv for f- at i,f+ at i.
A_LP = zeros(length_constraints,length_x_LP);
b_LP = zeros(length_constraints,1);
for i = 1:N_node
   for t = 1:N_task
       r_t = max(input_rate(:,t)); % rate of this task
       index_conserv_minus = (i-1)*N_task + t;
       index_conserv_plus = N_node * N_task + (i-1)*N_task + t;
       
% first set the conserv for f-
       b_LP(index_conserv_minus) = -1* input_rate(i,t);
       index_g = N_node^2 * 2 * N_task + (i-1)*N_task + t;
       A_LP(index_conserv_minus,index_g) = -1;
       for j = 1:N_node
           if Adj_cost(j,i) > 0 % in-going links (j,i)
               index_f_minus = (j-1)*N_node*N_task + (i-1)*N_task + t;
               A_LP(index_conserv_minus,index_f_minus) = 1;
           end
           if Adj_cost(i,j) > 0 % out-going links (i,j)
              index_f_minus = (i-1)*N_node*N_task + (j-1)*N_task + t;
              A_LP(index_conserv_minus,index_f_minus) = -1;
           end
       end

% then the conserv for f+
        if Task(t,1) ~= i
            index_g = N_node^2 * 2 * N_task + (i-1)*N_task + t;
            A_LP(index_conserv_plus,index_g) = a_m(Task(t,2)) +  + b_overhead / r_t;
            for j = 1:N_node
                if Adj_cost(j,i) > 0 % in-going links
                   index_f_plus = N_node^2*N_task + (j-1)*N_node*N_task + (i-1)*N_task + t;
                   A_LP(index_conserv_plus,index_f_plus) = 1;
                end
            end
        end
        for j = 1:N_node
           if Adj_cost(i,j) >0 % out-going links
               index_f_plus = N_node^2*N_task + (i-1)*N_node*N_task + (j-1)*N_task + t;
               A_LP(index_conserv_plus,index_f_plus) = -1;
           end
        end

   end
end

x_lb_LP = zeros(1,length_x_LP);
x_ub_LP = [UB*ones(1,N_node^2 * 2 * N_task) computation_cap*ones(1,N_node * N_task)];
[x_opt_LP, obj_opt_LP] = linprog(f_LP,A_LP,b_LP,[],[],x_lb_LP,x_ub_LP);
cost_opt_LP = 

%% compare
cost_opt_MILP
obj_opt_LP