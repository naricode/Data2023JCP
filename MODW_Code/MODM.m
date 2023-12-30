clc; clear; close all;
%% Problem Definition
CostFunction=@(x) MOTP02(x);      % Cost Function
nVar=3;             % Number of Decision Variables
VarSize=[1 nVar];   % Size of Decision Variables Matrix
VarMin=-4;          % Lower Bound of Variables
VarMax= 4;          % Upper Bound of Variables

% Number of Objective Functions
nObj=numel(CostFunction(unifrnd(VarMin,VarMax,VarSize)));

%% MODM-II Parameters
MaxIt=200;      % Maximum Number of Iterations
nPop=100;        % Population Size
nBabysitter= 3;         % Number of babysitters
nAlphaGroup=nPop-10*nBabysitter;         % Number of Alpha group
nScout=nAlphaGroup;         % Number of Scouts
L=round(0.6*nVar*nBabysitter); % Babysitter Exchange Parameter 
peep=2;             % Alpha female痴 vocalization 

global NFE
NFE=0;
%% Initialization
empty_mongoose.Position=[];
empty_mongoose.Cost=[];
empty_mongoose.Rank=[];
empty_mongoose.DominationSet=[];
empty_mongoose.DominatedCount=[];
empty_mongoose.CrowdingDistance=[];
pop=repmat(empty_mongoose,nAlphaGroup,1);

% Initialize Best Solution Ever Found
BestSol.Cost=inf;
tau=inf;
Iter=1;
sm=inf(nAlphaGroup,1);

% Create Initial Population
for i=1:nAlphaGroup
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    pop(i).Cost=CostFunction(pop(i).Position);
            if pop(i).Cost<=BestSol.Cost
             BestSol=pop(i);
            end
end
% Non-Dominated Sorting
[pop, F]=NonDominatedSorting(pop);

% Calculate Crowding Distance
pop=CalcCrowdingDistance(pop,F);

% Sort Population
[pop, F]=SortPopulation(pop);

% Store F1
%BestSol=[pop(F{1}).Cost];
% Abandonment Counter
C=zeros(nAlphaGroup,1);
CF=(1-Iter/MaxIt)^(2*Iter/MaxIt);

% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);

%% MODM Main Loop
tic
for it=1:MaxIt    
    %% Alpha group
     Ff=zeros(nAlphaGroup,1);
     MeanCost = mean(mean([pop.Cost]));
    for i=1:nAlphaGroup        
        % Calculate Fitness Values and Selection of Alpha
        W=mean([pop(i).Cost]);
        Ff(i) = exp(-W/MeanCost); % Convert Cost to Fitness
    end
        P=Ff/sum(Ff);
      % Foraging led by Alpha female
    for m=1:nAlphaGroup        
        % Select Alpha female
        i=RouletteWheelSelection(P);
        
        % Choose k randomly, not equal to Alpha
        K=[1:i-1 i+1:nAlphaGroup];
        k=K(randi([1 numel(K)]));
        
        % Define Vocalization Coeff.
        phi=(peep/2)*unifrnd(-1,+1,VarSize);     
        newpop.Position=[];
        newpop.Cost=[];
        newpop.Rank=[];
        newpop.DominationSet=[];
        newpop.DominatedCount=[];
        newpop.CrowdingDistance=[];

        % New Mongoose Position
        newpop.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
        
        % Evaluation
        newpop.Cost=CostFunction(newpop.Position);        
   
        % Comparision     
        if newpop.Cost<=pop(i).Cost
            pop(i)=newpop;
            % Non-Dominated Sorting
            [pop, F]=NonDominatedSorting(pop);

            % Calculate Crowding Distance
            pop=CalcCrowdingDistance(pop,F);

            % Sort Population
            [pop, F]=SortPopulation(pop);       
        else
            C(i)=C(i)+1;
        end    
    end     
    
    %% Scout group
    for i=1:nScout        
        % Choose k randomly, not equal to i
        K=[1:i-1 i+1:nAlphaGroup];
        k=K(randi([1 numel(K)]));
        
        % Define Vocalization Coeff.
        phi=(peep/2)*unifrnd(-1,+1,VarSize);        
        newpop.Position=[];
        newpop.Cost=[];
        newpop.Rank=[];
        newpop.DominationSet=[];
        newpop.DominatedCount=[];
        newpop.CrowdingDistance=[];        
        
        % New Mongoose Position
        newpop.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
        
        % Evaluation
        newpop.Cost=CostFunction(newpop.Position);              
   
        % Sleeping mould
        sm(i)=(mean(newpop.Cost)-mean(pop(i).Cost))/max(mean(newpop.Cost),mean(pop(i).Cost));
        
        % Comparision       
        if newpop.Cost<=pop(i).Cost
            pop(i)=newpop;
            % Non-Dominated Sorting
            [pop, F]=NonDominatedSorting(pop);

            % Calculate Crowding Distance
            pop=CalcCrowdingDistance(pop,F);

            % Sort Population
            [pop, F]=SortPopulation(pop);       
        else
            C(i)=C(i)+1;
        end       
    end   
    
    %% Babysitters
    for i=1:nBabysitter
        if C(i)>=L            
            pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
            pop(i).Cost=CostFunction(pop(i).Position);
            % Non-Dominated Sorting
            [pop, F]=NonDominatedSorting(pop);

            % Calculate Crowding Distance
            pop=CalcCrowdingDistance(pop,F);

            % Sort Population
            [pop, F]=SortPopulation(pop);
            C(i)=0;
        end
    end   
    
    % Update Best Solution Ever Found
      for i=1:nAlphaGroup
            if pop(i).Cost<=BestSol.Cost
             BestSol=pop(i);             
            % Non-Dominated Sorting
            [pop, F]=NonDominatedSorting(pop);

            % Calculate Crowding Distance
            pop=CalcCrowdingDistance(pop,F);

            % Sort Population
            [pop, F]=SortPopulation(pop);     
           end
      end    
        
   % Next Mongoose Position
   newtau=mean(sm);
   for i=1:nScout
        M=(pop(i).Position.*sm(i))/pop(i).Position;
        if newtau>tau
           newpop.Position=pop(i).Position-CF*phi*rand.*(pop(i).Position-M);
        else
           newpop.Position=pop(i).Position+CF*phi*rand.*(pop(i).Position-M);
        end
        tau=newtau;
   end       
  
    % Update Best Solution Ever Found
      for i=1:nAlphaGroup
            if pop(i).Cost<=BestSol.Cost
             BestSol=pop(i);
            % Non-Dominated Sorting
            [pop, F]=NonDominatedSorting(pop);

            % Calculate Crowding Distance
            pop=CalcCrowdingDistance(pop,F);

            % Sort Population
            [pop, F]=SortPopulation(pop);            
           end
      end
    
    % Store Best Cost Ever Found
    %BestCost(it)=BestSol.Cost;
    BEF=BestSol.Cost;
    BEP=BestSol.Position;
    BER=BestSol.Rank;
    BEC=BestSol.CrowdingDistance;
    F1=pop(F{1});  
        
    % Display Iteration Information
    disp(['Iteration = ' num2str(it) ' NFE = ' num2str(NFE) ' number of f1 = ' num2str(numel(F1))]);
    
    % Plot F1 Costs
    figure(1);
    PlotCosts(F1);  
end
Time=toc;
%name=char(datetime);
%name(name==':')='-';
%save(name)
Results.nPop=nPop;
Results.Iteration=MaxIt;
Results.Time=Time;
Results.NFE=NFE;
Results.F=numel(F1);
struct2table(Results)
%writetable(struct2table(Results),name)