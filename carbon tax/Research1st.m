% Clear the memory
clear all

%% Set experiment
reform = 1; % 0 -> without carbon taxation; 1  -> with carbon taxation;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demographics
J=81;                       % life-span(data)
JR=63;                      % age of retirement(data)
tR=J-JR+1;                  % length of retirement
tW=JR-1;                    % length of working life
n=0.004;                    % Population growth(data)
psi=1;                      % Survive prob(assumption)

% Preferences
beta=0.97;                  % discount factor（这个不知道怎么确定，但是大概这个数值也ok）(***)
theta1=2;                   % coefficient of relative risk aversion(Conesa et al. (2009))(***)
theta2=0.5;                 % Frisch Elasticity(***)
gamma=0.42;                 % weight on consumption(这个需要搜数据，看下家庭能源消费占家庭总消费的多少）(***)
Chi=1;                      % disutility of labor(assumption)



% Production
zeta=0.36;                  % Capital share(Data)(***)
phi=0.5;                    % Subsitution Elasticity(***)
alpha=0.36;                 % production elasticity of capital 
delta=0.06;                 % rate of depreciation（这个是否需要再去看看具体国家的情形是怎样的）
A=1;                        % Productivity(Normalize)
Pe=0.0025;                  % Energy price(data)(***)


% Measure of each generation
mass=ones(J,1);
for ik0=2:J
    mass(ik0)=mass(ik0-1)/(1+n);
end

% Normalized measure of each generation (sums up to 1)
mass=mass/sum(mass);

% Age-efficiency profile(***)(这个需要进一步讨论如何做）
e = load('ef.txt');
plot(e)

%  Capital grid
maxkap = 14;                               % maximum value of capital grid  
minkap = 0.01;                             % minimum value of capital grid
nk=800;                                    % number of grid points
inckap=(maxkap-minkap)/(nk-1);             % distance between points
aux=1:nk;
kap= minkap+inckap*(aux-1);                % capital grid

neg=-1e10;                                 % very small number

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Carbon taxation and initial guesses(Make a guess on equilibrium quantity of assets and effective labor)
if reform == 1
    tauC=50;
    K0=3.9875;           
    N0=0.3800;
    E0=8000; 
else
    tauC=0;
    K0= 2.9734;
	N0=0.3683;
    E0=9000;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop over capital and labor(OLG loop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tolk=1e-3;              % Numerical tolerance for capital（###）
tollab=1e-3;            % Numerical tolerance for labor
nq=10;                  % Max number of iterations

q=0;                    % Counter for iterations
K1=K0+10;
N1=N0+10;

fprintf('\nComputing equilibrium price... \n');
while q<nq && (abs(K1-K0)>tolk || abs(N1-N0)>tollab)
   
    q=q+1;
    
    fprintf('\nIteration %g out of %g \n',q,nq);
    
     % Prices
    r0 = zeta*((((K0^zeta)*(N0^(1-zeta)))^((phi-1)/phi)+(E0)^((phi-1)/phi))^(1/(phi-1)))*((K0^zeta)*(N0^(1-zeta)))^(-1/phi)+(E0)^(1-zeta)-delta;
    w0 = (1-zeta)*(((((K0^zeta)*(N0^(1-zeta)))^((phi-1)/phi)+(E0)^((phi-1)/phi))^(1/phi-1)))*((K0^zeta)*(N0^(1-zeta)))^(-1/phi)*(K0/N0)^zeta;
   
    % Carbon taxation revenue
    cr = sum((tauC*ec)*mass);  %这个公式不确定
    

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BACKWARD INDUCTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initialization
    v=zeros(nk,J);                          % value function of agents(two dimensions:one is for capital number, one is for age)
    kapopt=ones(nk,J);                      % optimal savings of agents
    % (store INDEX of k' in capital grid, not k' itself!)
    labopt=ones(nk,tW);                     % optimal labor supply   
    
        
    % Retired households
    
    % Last period utility
    cons=(1+r0)*kap-(Pe+tauC)*ec;                    % last period consumption (vector!)(***)（总是觉得这里有点奇怪）
    util=((cons.^gamma*ec.^(1-gamma))^(1-theta1))/(1-theta1);       % last period utility (vector!)(.:向量乘积）(l=0,since is the last period)
    v(:,J)=util;                            % last period indirect utility (vector!)
    
    for j=J-1:-1:tW+1 % age(-1:means backwards)(tW+1=JR)
        for ik0=1:nk        % assets today（***）
   %   We need to solve for optimal savings decision at age j=J-1 for each(line129:ik0=1:nk, loop in line) value of current  asset holdings        
            % Initialize right-hand side of Bellman equation
            vmin=neg;
            ik1=0;
            
            %agent with current asset aJ-1=0(1st guy from asset grid)
            % Idea:take a1J(first element in asset grid), then take the
            % second element,until the 800th element. Picl aj that
            % maximizes RHS of the Bellman equation and store it.
            
            % Loop over all k's in the capital grid to find the value,
            % which gives max of the right-hand side of Bellman equation
            
            while ik1<nk  	% assets tomorrow
                ik1=ik1+1; 
                kap0=kap(ik0); % current asset holdings
                kap1=kap(ik1); % future asset holdings
                
                % Instantaneous utility
                cons=(1+r0)*kap0-kap1-(Pe+tauC)*ec;
                
                if cons<=0
                    util=neg;
                else
                    util=((cons.^gamma*ec.^(1-gamma))^(1-theta1))/(1-theta1);
                end
                
                % Right-hand side of Bellman equation
                v0=util + beta*v(ik1,j+1);
                
                % Store indirect utility and optimal saving
                if v0>vmin
                    v(ik0,j)=v0;
                    kapopt(ik0,j)=ik1;
                    vmin=v0;
                end
            end
        end
    end
    
    % Working households
    for j=tW:-1:1           % age
        for ik0=1:nk          % assets today
            
            % Initialize right-hand side of Bellman equation
            vmin=neg;
            ik1=0;
            
            % Loop over all k's in the capital grid to find the value,
            % which gives max of the right-hand side of Bellman equation
            
            while ik1 < nk  	% assets tomorrow
                ik1=ik1+1;
                
                kap0=kap(ik0); % current asset holdings
                kap1=kap(ik1); % future asset holdings
                
                % Optimal labor supply
                % CHECK THIS!
                lab=(gamma*(1-tau)*e(j)*w0-(1-gamma)*((1+r0)*kap0-kap1))/((1-tau)*e(j)*w0);%需要修改，energy与consumption的关系不明晰
                
                % Check feasibility of labor supply
                if lab>1
                    lab=1;
                elseif lab<0
                    lab=0;
                end
    
  % Instantaneous utility
                cons=(1+r0)*kap0+(1-tau)*w0*e(j)*lab-kap1-(Pe+tauC)*ec;
                
                if cons<=0
                    util=neg;
                else
                    util=(((cons^gamma)*((ec)^(1-gamma)))^(1-theta1))/(1-theta1)-(h^(1+(1/theta2)))/(1+(1/theta2));
                end
                
                % Right-hand side of Bellman equation   
                v0=util + beta*v(ik1,j+1);
                
                % Store indirect utility, optimal saving and labor
                if v0>vmin
                    v(ik0,j)=v0;
                    kapopt(ik0,j)=ik1;
                    labopt(ik0,j)=lab;
                    vmin=v0;
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Aggregate capital stock and employment                                  % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initializations
    kgen=zeros(J,1);        % capital supply k(j+1) for each generation j
    labgen=zeros(tW,1);     % labor supply l(j) for each generation j

    % Use decision rules to iteratively fill in kgen and labgen
    ik0 = 1;                % starting capital of j = 1, kap(ik0) = 0(***)
    for j=1:J               % iterations over cohort
        % capital decision kp(k)(to find the optimal one)
        ik1 = kapopt(ik0, j);
        kgen(j) = kap(ik1);
        
        % labor decision l(k)
        if j<=tW
            labgen(j) = labopt(ik0, j);
        end

        % update k = kp(***)
        ik0 = ik1;
    end

    K1 = kgen' * mass;              % dot product of vectors
    L1 = (labgen .* e)' * mass(1:tW);      % dot product of vectors

    % Update the guess on capital and labor    
    K0=0.9*K0+0.1*K1;
    N0=0.9*N0+0.1*N1;

    % Display results
    disp('  capital     labor   carbon taxation');
    disp([K0, N0, ]);
    disp('deviation-capital deviation-labor       ');
    disp([abs(K1-K0),  abs(N1-N0)]);
end 


% Average hours worked
h = labgen'* mass(1:tW)/sum(mass(1:tW));

% Gini disposable income
income = zeros(J,1);
income(1:tW) = (1-tau)*labgen .* e*w0 + r0*kgen(1:tW)-(Pe+tauC)*ec;% ec不知道怎么体现(***)
income(tW+1:end) = r0*kgen(tW+1:end);
gini_index = gini(mass,income);


% Display equilibrium results
% Output
Y=((K0^zeta*N0^(1-zeta))^((phi-1)/phi)+(Ep)^((phi-1)/phi))^(phi/(phi-1));

disp('      K0         N0       w         r         carbon taxation     Y     K/Y    h    Gini');
disp([K0, L0, w0, r0, b, Y, K0/Y, h, gini_index]);

%% Comparisons across steady-states
if reform ==0
    save('ss.mat');
else
	save('no_ss.mat');
end



