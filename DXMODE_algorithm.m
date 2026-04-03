% DXMODE_algorithm.m
% Official MATLAB implementation of the DXMODE algorithm
%
% Paper:
% Reda, M., Onsy, A., Haikal, A.Y., and Ghanbari, A.
% DXMODE: A dynamic explorative multi-operator differential evolution
% algorithm for engineering optimization problems.
% Information Sciences, 717 (2025), 122271.
% DOI: 10.1016/j.ins.2025.122271
%
% Inputs:
%   CECyear : benchmark year (e.g., 2020 or 2022)
%   fNo     : benchmark function number
%   nd      : problem dimension
%   lb      : lower bound (scalar or 1-by-nd vector)
%   ub      : upper bound (scalar or 1-by-nd vector)
%
% Outputs:
%   goalReached : logical flag indicating whether tolerance was reached
%   GlobalBest  : best solution structure with fields Position and Cost
%   countFE     : total number of function evaluations

function [goalReached, GlobalBest, countFE] = DXMODE_algorithm(CECyear, fNo, nd, lb, ub)

    %% Begin intialiazation
    %% CEC function paramters

    % Default values if no inputs are passed
    if nargin < 1 || isempty(CECyear)
        CECyear = 2020;
    end

    if nargin < 2 || isempty(fNo)
        fNo = 3;
    end

    if nargin < 3 || isempty(nd)
        nd = 20;
    end

    if nargin < 4 || isempty(lb)
        lb = -100;
    end

    if nargin < 5 || isempty(ub)
        ub = 100;
    end

    %% Cost Function
    if CECyear == 2020
        costFunction = @(sol) cost_cec2020(sol, fNo);
    elseif CECyear == 2022
        costFunction = @(sol) cost_cec2022(sol, fNo);
    else
        error('Unsupported CEC year. Use 2020 or 2022.');
    end


    %% Algorithm paramters
    AlgName = 'DXMODE Algorithm' ;
    n_opr = 4 ; % number of mutation operators 
    
    nPopRate = 12 ; % initial population is set to nPopRate * nd 
    MinPopSize = 6 ; % min pop size (linear reduction for population)

    arch_rate = 2.6; % 2.6 (A) rate of the archive control the size of it
    mem_size_scale = 6 ; % memory size (H) = 20 * nd, where 20 is the memeory scale
    
    Cr_init = 0.2; % Initial CR
    F_init = 0.2; % initial F

    PX_max = 0.4; %0.4  %0.25
    PX_min = 0.05; % 0.05

    maxfe_ls_ratio = 0.85; % ration at which the local search starts (after 85% of maxfes)
    max_ls_iterations_rate = 20.0000e-003; % CFE_ls max iterations for ls (ration of maxfe)
    prob_ls=0.1;  % porabibility of quadratic programming optimization using fmincon
    prob_ls_min = 0.01 ; % 0.0001 in the paper 0.01 in code set this min probability if the LS not improving the global vest
    ls_optimizar_name = 'sqp' ; % sequantial quadratic programming 'sqp' or  'active-set' or 'trust-region-reflective' or 'interior-point'

    n_opr_explor = 3 ; % number of explration operators (methods)
    sigma = 0.1; % Standard deviation of the Gaussian distribution, adjust as needed
    levy_exploration_prob = 0.25 ; % explration probabilty of levy random walk (Pa)
    Chaotic_max_loop_count = 10; % max number of iteration for performing chaotic map update on a solution                   
    Chaotic_mapTypes = {'logistic', 'sine', 'tent'};  % List of available chaotic map types
   
    % intial pop size 
    nPop = nPopRate * nd  ; % Initial Population Size
    InitPop = nPop ;
    % prob. of each DE operator
    probDE1= 1./n_opr .* ones(1,n_opr);
    prob_explor = 1./n_opr_explor .* ones(1,n_opr_explor);

    %% Initialize Archive Data 
    archive.NP = arch_rate * nPop; % the maximum size of the archive
    cand = create_empty_individual(); 
    cand.Position = zeros(1, nd); % empty solution of nd dim
    cand.Cost = inf(1, 1); % empty cost value
    archPop = repmat(cand, 1, 1);  % popualtion of one individual
    archive.pop = archPop; 

    %% Initialize Adaptive Archive for CR and F
    hist_pos=1;
    memory_size=mem_size_scale*nd;
    archive_f = ones(1,memory_size).* F_init;
    archive_Cr = ones(1,memory_size).* Cr_init;

    %% Stopping criteria
    tol = 10^-8;
    if CECyear == 2020 || CECyear == 2022
        if nd == 10
            maxfe = 200000 ; 
        elseif nd == 20 
            maxfe = 500000;
        else
            fprintf('Dimensions must be 10 or 20 \n');
            return
        end
    elseif CECyear == 2017
        maxfe = nd * 10000;
    end

    %% Display iteration prompt
    print_flag = true;

    %%  Global variable to count number of function evaluations
    global countFE;
    countFE = 0 ;

    %% Initialize iteration counter 
    N_iter = 0;

    %% Goal reached flag
    goalReached = false; 

    %% Initialize Global best
    GlobalBest.Position = [];
    GlobalBest.Cost = Inf ; 

    %% Set the seed for random number generator
    rng('default');  % Resets to the default settings
    rng('shuffle'); % set it to shuffle
    
    %% Initialize population and update the global best
    population = repmat(create_empty_individual(), nPop, 1);

    % Generate initial positions using Latin Hypercube Sampling
    LB = repmat(lb, 1, nd);
    UB = repmat(ub, 1, nd);
    LHS_samples = lhsdesign(nPop, nd); % Generates an nPop x CEC.dim matrix of samples
    % Scale samples to the problem's bounds
    for d = 1:nd
        LHS_samples(:, d) = LB(d) + (UB(d) - LB(d)) .* LHS_samples(:, d);
    end
    
    for i = 1:nPop
        % Initialize Position with scaled LHS samples
        if rand <=1
            population(i).Position = LHS_samples(i, :);
        else
            population(i).Position = lb + (ub - lb) .* rand(1, nd);
        end
        
        % Evaluation of the cost
        population(i).Cost = costFunction(population(i).Position);

        % Update Global Best
        if population(i).Cost < GlobalBest.Cost
            GlobalBest = population(i);
        end
    end

    % Initial Error (Global Best Cost)
    InitError = GlobalBest.Cost ; % used for population size reduction
    % create a random permuatation of the popualtion
    pop_old = population(randperm(nPop),:);

   %% begin algorithm loop 
    while (GlobalBest.Cost > tol)  && (countFE <= maxfe)
        %% update the generation
        N_iter=N_iter+1; 

        %% Update popuation size
        EIR = 1 - ((GlobalBest.Cost - tol) / (InitError - tol));
        FERate = (countFE / maxfe) ; 
        TR = 0.4 * EIR + 0.6 * FERate; 
        newPopSize= round(MinPopSize + ((InitPop - MinPopSize) * (1 - ( TR).^1)));

        %% Update the popuation according to the new popuation size
        nPop = numel(population); % current popsize 
        if nPop > newPopSize
            % Calculate the number of individuals to remove
            reduction_ind_num = nPop - newPopSize;
            if nPop - reduction_ind_num < MinPopSize
                reduction_ind_num = nPop - MinPopSize;
            end
        
            % Remove the worst individuals
            for r = 1 : reduction_ind_num
                % Sort population based on Cost
                [~, sortedIdx] = sort([population.Cost], 'descend');
                % Remove the worst individual
                population(sortedIdx(1)) = []; % it removed from the original popualtion
            end

            % update the current popsize
            nPop = numel(population);
        
            %% Update archive size based on the new population size
            archive.NP = round(arch_rate * nPop);

            % If archive size exceeds its limit, randomly remove some individuals
            current_archive_NP = numel(archive.pop);
            if current_archive_NP > archive.NP
                rndpos = randperm(current_archive_NP);
                rndpos = rndpos(1 : archive.NP);
                archive.pop = archive.pop(rndpos);
            end
        end

        %% Initialize the archive of the CR and F
        mem_rand_index = ceil(memory_size * rand(nPop, 1));
        mu_sf = archive_f(mem_rand_index);
        mu_cr = archive_Cr(mem_rand_index);
        
        %%  generate CR   
        cr = normrnd(mu_cr, 0.1);
        term_pos = find(mu_cr == -1);
        cr(term_pos) = 0;
        cr = min(cr, 1);
        cr = max(cr, 0);
        % sort the cr
        [cr,~]=sort(cr);

        %% for generating scaling factor
        F = mu_sf + 0.1 * tan(pi * (rand(1,nPop) - 0.5));
        pos = find(F <= 0);
        
        while ~ isempty(pos)
            F(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand(1,length(pos)) - 0.5));
            pos = find(F <= 0);
        end
        
        F = min(F, 1);
        F=F';

         %% Sort the popuation     
        Costs = [population.Cost]; % original costs of the original popuation
        [Costs, SortOrder] = sort(Costs);
        population = population(SortOrder);  

        %% **** Mutation Phase ****
        % combine the popuation with the archive population  
        popAll = [population; archive.pop];  
        
        %% generate mutation operator probablities for each individual in the population
        % Randomly decide the mutation strategy for each individual
        bb = rand(nPop, 1);
        
        % Retrieve probabilities for each strategy
        probiter = probDE1(1, :);
        l2 = sum(probDE1(1:2));
        l3 = sum(probDE1(1:3)); 
        % Determine which strategy to apply for each individual
        op_1 = bb <= probiter(1) * ones(nPop, 1);
        op_2 = (bb > probiter(1)) & (bb <= l2);
        op_3 = (bb > l2) & (bb <= l3);
        op_4 = (bb > l3) & (bb <= 1);

        %% generate random integer numbers
        r0 = 1 : nPop;
        [r1, r2,r3] = gnR1R2(nPop, size(popAll, 1), r0);

        %% Choose top individuals (at least one) for DE operator 1 and 2
        pNP12 = max(round(0.25 * nPop), 1); % At least one or 25% of the population size
        randindex = ceil(rand(1, nPop) .* pNP12); % Select indices from the best subset
        randindex = max(1, randindex); % Ensuring indices are valid (not less than 1)
        phix12 = population(randindex, :);
        
        %%  Choose top individuals (at least two) for DE operator 3
        pNP3 = max(round(0.5 * nPop), 2); %% choose at least two best solutions
        randindex = ceil(rand(1, nPop) .* pNP3); %% select from [1, 2, 3, ..., pNP]
        randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
        phix3 = population(randindex, :);

        %% Initialize mutation vector
        cand = create_empty_individual();
        cand.Position = zeros(1, nd);
        newPop = repmat(cand, nPop, 1);

        for i = 1:nPop
            % apply mutation rule
            x_curr = population(i).Position;
            x_r1 = population(r1(i)).Position;
            x_r3 = population(r3(i)).Position;
            xx_r2 = popAll(r2(i)).Position;
            x_phi12 = phix12(i).Position;
            x_phi3 = phix3(i).Position;

            if op_1(i)
                % Strategy 1: DE/current-to-pbest/1/Arch
                newPop(i).Position = x_curr + F(i) * (x_phi12 - x_curr + x_r1 - xx_r2);
            elseif op_2(i)
                % Strategy 2: DE/current-to-pbest/1
                newPop(i).Position = x_curr + F(i) * (x_phi12 - x_curr + x_r1 - x_r3);
            elseif op_3(i)
                % Strategy 3: DE/weighted-rand-to-qbest/1
                newPop(i).Position = F(i) * x_r1 + F(i) * (x_phi3 - x_r3);
            elseif op_4(i)
                % Strategy 4: DE/rand/1
                idxs=randperm(nPop);    
                idxs(idxs==i)=[];  %exclude the current position  
                a = population(idxs(1)).Position;
                b = population(idxs(2)).Position;
                c = population(idxs(3)).Position; 
                newPop(i).Position = a + F(i).* (b - c);
            end
             % Applying boundary check for each individual
             newPop(i).Position = han_boun_individual(newPop(i).Position, ub, lb, x_curr);
        end

        %% *** Crossover ***
        % Initialize ui as an array of individuals
        cand = create_empty_individual();
        cand.Position = zeros(1, nd);
        newPop2 = repmat(cand, nPop, 1);

        for i = 1:nPop

            % Binomial Crossover
            mask = rand(1, nd) > cr(i);
            jrand = floor(rand * nd) + 1; % Ensure at least one dimension is inherited from vi
            mask(jrand) = false;
            newPop2(i).Position = population(i).Position; % Start with parent position from origial population
            newPop2(i).Position(~mask) = newPop(i).Position(~mask); % Inherit from newPop where mask is false

            %% evaluate the new cost
            newPop2(i).Cost = costFunction(newPop2(i).Position);   
        end

        %% *** Update the archives ****
        %% get the I label for the improved individuals
        newCosts = [newPop2.Cost] ;
        I = (newCosts < Costs ); % Logical index of improved solutions

        %% update the archive with the old bad solutions in the population
        archive = updateArchive(archive, population(I == 1));

        %% update probDE1 (operators probabilities). of each DE 
        diff2 = max(0, (Costs - newCosts))./abs(Costs + eps); % Improvement metric, adding eps for stability

        % Calculate performance scores for this iteration
        count_S = zeros(1, n_opr);
        count_S(1)=max(0,mean(diff2(op_1==1)));
        count_S(2)=max(0,mean(diff2(op_2==1)));
        count_S(3)=max(0,mean(diff2(op_3==1)));
        count_S(4)=max(0,mean(diff2(op_4==1)));

        % Check if there is any significant improvement across all operators
        if all(count_S <= eps) % If no significant improvement, reset to equal probabilities
            probDE1 = ones(1, n_opr) / n_opr;
        else
            probDE1 = max(0.1, min(0.9, probDE1));
        end

        %% calc. imprv. for Cr and F archives (as in Ref [9])
        goodCR = cr(I == 1);
        goodF = F(I == 1);
        diff = abs(Costs - newCosts);
        if size(goodF,1)==1
            goodF=goodF';
        end
        if size(goodCR,1)==1
            goodCR=goodCR';
        end
        num_success_params = numel(goodCR);
        if num_success_params > 0
            weightsDE = diff(I == 1)./ sum(diff(I == 1));
            %% for updating the memory of scaling factor
            archive_f(hist_pos) = (weightsDE * (goodF .^ 2))./ (weightsDE * goodF);
            
            %% for updating the memory of crossover rate
            if max(goodCR) == 0 || archive_Cr(hist_pos)  == -1
                archive_Cr(hist_pos)  = -1;
            else
                archive_Cr(hist_pos) = (weightsDE * (goodCR .^ 2)) / (weightsDE * goodCR);
            end
            
            hist_pos= hist_pos+1;
            if hist_pos > memory_size;  hist_pos = 1; end
        else
            archive_Cr(hist_pos)=0.5;
            archive_f(hist_pos)=0.5;
        end

        %% update population with the good solution and update the oldPop with the bad old solutions
        pop_old(I == 1) = population(I == 1); %save the bad individual in popualtion in the old popuation
        population(I == 1) = newPop2(I == 1); % relace the bad individuals in popuation with the better individuals in newPop2

        %% sort the population and old population
        [~, sortedIndices] = sort([population.Cost]); % the new updated merged costs
        population = population(sortedIndices);
        pop_old = pop_old(sortedIndices);

        %% update global best
        localBest = population(1); % the sorted population, the top solution is the best cost solution
        if localBest.Cost < GlobalBest.Cost  
            GlobalBest = localBest;
        end

        %% *** random exploration ***
        prob_explore = PX_max - (PX_max - PX_min) * (countFE/maxfe);

        if rand < prob_explore
           oldCosts = [population.Cost] ;
    
            %% Initialize new population 
            cand = create_empty_individual();
            cand.Position = zeros(1, nd);
            newPop3 = repmat(cand, nPop, 1);
    
            %% select exporation operator probablities for each individual in the population
            % Randomly decide the mutation strategy for each individual
            bb = rand(nPop, 1);
          
            % Retrieve probabilities for each strategy
            probiter = prob_explor(1, :);
            l2 = sum(prob_explor(1:2));
    
            % Determine which strategy to apply for each individual
            op_e1 = bb <= probiter(1) * ones(nPop, 1);
            op_e2 = (bb > probiter(1)) & (bb <= l2);
            op_e3 = (bb > l2) & (bb <= 1);

            cc = 0 ; 
            for i = 1:nPop
                
               newSol = create_empty_individual();
               mapId = randi(length(Chaotic_mapTypes)) ; % {'logistic', 'sine', 'tent'}

                if op_e1(i)  
                    % Strategy 1: levey random walk
                    newSol.Position = BSRW_exp(i,population, levy_exploration_prob  );
                elseif op_e2(i)
                    % Strategy 2: Chaotic maps 
                    newSol.Position = MNCE_chaotic(i,population, Chaotic_mapTypes, Chaotic_max_loop_count, mapId);
                   
                elseif op_e3(i)
                    % Strategy 3: Gaussian Mutation
                    newSol.Position = population(i).Position + sigma * randn(size(population(i).Position));
                end
                oldSol = newSol.Position ;  

               % handle boundaries
               newSol.Position = han_boun_individual(newSol.Position, ub, lb, population(i).Position);
    
                if isequal(newSol.Position, oldSol)
                    cc = cc + 1 ;    
                end
                newPop3(i).Position = newSol.Position;

                %% evaluate the new cost
                newPop3(i).Cost = costFunction(newPop3(i).Position);
            end
           newCosts = [newPop3.Cost] ;
    
           I2 = (newCosts < oldCosts ); % Logical index of improved solutions

            %% update the archive with the old bad solutions in the population
            archive = updateArchive(archive, population(I2 == 1)) ; 
    
          %% update exploration  (operators probabilities).
            diff3 = max(0, (oldCosts - newCosts))./abs(oldCosts + eps); % Improvement metric, adding eps for stability
  
            % Calculate performance scores for this iteration
            count_ex_S = zeros(1, n_opr_explor);
            count_ex_S(1)=max(0,mean(diff3(op_e1==1)));
            count_ex_S(2)=max(0,mean(diff3(op_e2==1)));
            count_ex_S(3)=max(0,mean(diff3(op_e3==1)));
    
            % Check if there is any significant improvement across all operators
            if all(count_ex_S <= eps) % If no significant improvement, reset to equal probabilities
                prob_explor = ones(1, n_opr_explor) / n_opr_explor;
            else
                prob_explor = max(0.1, min(0.9, prob_explor));
            end
            
            %% update population with the good solution and update the oldPop with the bad old solutions
            pop_old(I2 == 1) = population(I2 == 1); %save the bad individual in popualtion in the old popuation
            population(I2 == 1) = newPop3(I2 == 1); % relace the bad individuals in popuation with the better individuals in newPop2
    
            %% sort the population and old population
            [~, sortedIndices] = sort([population.Cost]); % the new updated merged costs
            population = population(sortedIndices);
            pop_old = pop_old(sortedIndices);
    
            %% update global best
            localBest = population(1); % the sorted population, the top solution is the best cost solution
            if localBest.Cost < GlobalBest.Cost  
                GlobalBest = localBest;
            end
     
        end % end of exploration phase

        %% *** Local Search ***
        if countFE > maxfe_ls_ratio * maxfe && countFE<maxfe
            if rand<prob_ls
                [GlobalBest, succ] = LS2(GlobalBest, ub, lb, nd, maxfe, ls_optimizar_name, max_ls_iterations_rate, costFunction );
                if succ==1 %% if LS2 was successful
                    % replace the worst solution in the popualtion, with the
                    % new global best, because the new global best is not in  the orgignial population
                    population(nPop) = GlobalBest;
    
                    % sort the population and old population
                    [~, sortedIndices] = sort([population.Cost]); % the new updated merged costs
                    population = population(sortedIndices);
                    pop_old = pop_old(sortedIndices);
    
                    prob_ls = prob_ls;
                else
                    prob_ls=prob_ls_min; %% set p_LS to a small value it  LS was not successful
                end    
            end
        end

        %% check if maxfes is exceeded 
        if countFE > maxfe 
            break;
        end

        %% print the iteration number
        if print_flag            
            fprintf('%s | CEC%d_F%d_D%d | Iteration %d |FEs %d | Error %d\n', AlgName,  CECyear, fNo , nd, N_iter, countFE,GlobalBest.Cost);
        end

         %check tolerance /error
        if (GlobalBest.Cost <= tol)
            GlobalBest.Cost  = 0 ;
            disp('tol reached');
            goalReached = true ; 
            break ;  % not needed, becuase it will exit in the next while loop check
        end

    end
end

function individual = create_empty_individual()
    individual.Position = [];
    individual.Cost = Inf;
end

function archive = updateArchive(archive, newPop)
    % Update the archive with input solutions
    % Input:
    %   archive - The existing archive with fields 'NP' and 'pop'
    %   newPop - Array of new individuals to add to the archive

    if archive.NP == 0, return; end

    % Combine existing archive population with new population
    combinedPop = [archive.pop; newPop];

    % Extract positions for comparison
    positions = arrayfun(@(x) x.Position, combinedPop, 'UniformOutput', false);

    % Remove duplicates in positions
    [~, uniqueIdx] = unique(cell2mat(positions), 'rows');
    combinedPop = combinedPop(uniqueIdx);

    % Randomly remove solutions if necessary to maintain archive size
    if numel(combinedPop) > archive.NP
        rndIdx = randperm(numel(combinedPop), archive.NP);
        archive.pop = combinedPop(rndIdx);
    else
        archive.pop = combinedPop;
    end
end



function [GlobalBest, succ] = LS2(GlobalBest, ub, lb, nd, Max_FES, ls_optimizar_name, max_ls_iterations_rate, costFunction )
    global countFE; 
    
    % set the upper and lower bounds in the size of (nd *1)
    ub = repmat(ub, nd, 1);
    lb = repmat(lb, nd , 1);
    
    % set the maxfes for the fmincon algrotihm (small ration of the original maxfes)
    LS_FE = min(ceil(max_ls_iterations_rate * Max_FES), (Max_FES - countFE)); 
   
    % set the option for the optimizaer
    options = optimset('Display', 'off', 'algorithm', ls_optimizar_name, 'UseParallel', 'never', 'MaxFunEvals', LS_FE);

    % Adapted to use the new cost function format
    costFuncWrapper = @(X) costFunction(X');  % Note the transpose of X
    
    % run the algorithm
    [Xsqp, newCost, ~, details] = fmincon(costFuncWrapper, GlobalBest.Position', [], [], [], [], lb, ub, [], options);

    % Check if there is an improvement in the fitness value and update P_{ls}
    if (GlobalBest.Cost - newCost) > 0
        succ = 1;
        GlobalBest.Position = Xsqp';
        GlobalBest.Cost = newCost ;
    else
        succ = 0;
    end

end

function [r1, r2,r3] = gnR1R2(NP1, NP2, r0)
  
    NP0 = length(r0);
    r1 = floor(rand(1, NP0) * NP1) + 1;
    
    for i = 1 : 99999999
        pos = (r1 == r0);
        if sum(pos) == 0
            break;
        else % regenerate r1 if it is equal to r0
            r1(pos) = floor(rand(1, sum(pos)) * NP1) + 1;
        end
        if i > 1000, % this has never happened so far
            error('Can not genrate r1 in 1000 iterations');
        end
    end
    
    r2 = floor(rand(1, NP0) * NP2) + 1;
    %for i = 1 : inf
    for i = 1 : 99999999
        pos = ((r2 == r1) | (r2 == r0));
        if sum(pos)==0
            break;
        else % regenerate r2 if it is equal to r0 or r1
            r2(pos) = floor(rand(1, sum(pos)) * NP2) + 1;
        end
        if i > 1000, % this has never happened so far
            error('Can not genrate r2 in 1000 iterations');
        end
    end
    
    r3= floor(rand(1, NP0) * NP1) + 1;
    %for i = 1 : inf
    for i = 1 : 99999999
        pos = ((r3 == r0) | (r3 == r1) | (r3==r2));
        if sum(pos)==0
            break;
        else % regenerate r2 if it is equal to r0 or r1
             r3(pos) = floor(rand(1, sum(pos)) * NP1) + 1;
        end
        if i > 1000, % this has never happened so far
            error('Can not genrate r2 in 1000 iterations');
        end
    end
end

function x = han_boun_individual(x, ub, lb, x2)
    if isscalar(ub)
        ub = repmat(ub, 1, numel(x2));
    end
    if isscalar(lb)
        lb = repmat(lb, 1, numel(x2));
    end

    x_L = lb;
    pos = x < x_L;
    x(pos) = (x2(pos) + x_L(pos)) / 2;

    x_U = ub;
    pos = x > x_U;
    x(pos) = (x2(pos) + x_U(pos)) / 2;
end

function newPos = chaoticMap(pos, mapType)
    % Chaotic Map Function
    % Inputs:
    %   pos - Current position
    %   mapType - Type of chaotic map ('logistic', 'sine', 'tent')
    % Output:
    %   newPos - New position after applying chaotic map

    if nargin < 2
        mapType = 'logistic'; % Default map
    end

    switch mapType
        case 'logistic'
            % Logistic Map
            r = rand * (3.57- 0.1) + 4; % Random r in [3.57, 4]
            newPos = r .* pos .* (1 - pos);
        case 'sine'
            % Sine Map
            r = rand * (1 - 0.9) + 0.9; % Random r in [0.9, 1]
            newPos = r .* sin(pi .* pos);
        case 'tent'
            % Tent Map
            r = rand * (5 - 1.5) + 1.5; % Random r in [1.5, 2]
            if pos < 0.5
                newPos = r .* pos;
            else
                newPos = r .* (1 - pos);
            end
        otherwise
            error('Unknown map type');
    end
end

function newPos = MNCE_chaotic(i,population, Chaotic_mapTypes, Chaotic_max_loop_count, mapId)
        % Strategy 2 : Chaotic maps 
        % Randomly select a map type
        selectedMapType = Chaotic_mapTypes{mapId};
                      
        % Randomly select an integer value between 1 and nLoop  for repeatintion fo the map
        ll = randi(Chaotic_max_loop_count);
        newPos = chaoticMap(population(i).Position, selectedMapType);
        for l = 1: ll
           newPos = chaoticMap(newPos, selectedMapType);
        end

end


function newPos = BSRW_exp(i,population, levy_exploration_prob  )
        % Generate a random matrix K for selective modification
        nPop = numel(population);
        K = rand(size(population(i).Position)) > levy_exploration_prob;

        % Randomly select two different nests for step size calculation
        idx = randperm(nPop, 2);
        stepsize = rand * (population(idx(1)).Position - population(idx(2)).Position);

        % Update the nest position selectively
        newPos = population(i).Position + stepsize .* K;
end





