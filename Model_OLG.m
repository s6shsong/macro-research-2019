% Clear the memory
clear all

%% Set experiment
reform = 1; % 0 -> without carbon taxation; 1  -> with carbon taxation;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demographics
J=81;                       % life-span(data)
JR=65;                      % age of retirement(data)
tR=J-JR+1;                  % length of retirement
tW=JR-1;                    % length of working life
n=0.0047;                    % Population growth(data)
psi=1;                      % Survive prob(assumption)(simplify)

% Preferences
beta=0.97;                  % discount factor
theta1=2;                   % coefficient of relative risk aversion(Conesa et al. (2009))
theta2=0.5;                 % Frisch Elasticity
gamma=0.42;                 % weight on consumption
Chi=1;                      % disutility of labor(assumption)

% Production
zeta=0.36;                  % Capital share(Data)
phi=0.5;                    % Subsitution Elasticity
delta=0.083;                 % rate of depreciation
A=1;                        % Productivity
Pe=0.0025;                  % Energy price(Normalize)

% Government parameters
tauL=0.439;                  % Labor tax rate(Data)
G = 0.443;                   % government spending(normalize)


% Measure of each generation
mass=ones(J,1);
for ik0=2:J
    mass(ik0)=mass(ik0-1)/(1+n);
end

% Normalized measure of each generation (sums up to 1)
mass=mass/sum(mass);

% Age-efficiency profile(There will be two profiles if we divide them in to el and eh)
e = load('e.txt');
plot(e)

%  Capital grid
maxkap = 14;                               % maximum value of capital grid  
minkap = 0.01;                             % minimum value of capital grid
nk=800;                                    % number of grid points
inckap=(maxkap-minkap)/(nk-1);             % distance between points
aux=1:nk;
kap= minkap+inckap*(aux-1);                % capital grid
neg=-1e10;                                % very small number


% Energy consumption grid

maxene = 20;                              % maximum value of energy consumption grid
minene = 0.5;                             % minimum value of capital grid
nk=800;                                   % number of grid points
incene=(maxene-minene)/(nk-1);            % number of grid points
aux=1:nk;
ene= minene+incene*(aux-1);               % energy consumption grid


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Carbon taxation and initial guesses(Make a guess on equilibrium quantity of assets and effective labor)
if reform == 1
    tauC=50;
    K0=3.9875;           
    E0=3.543; 
    N0=0.3800;
else
    tauC=0;
    K0= 2.9734;
    E0=4.387;
	N0=0.3683;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop over capital，energy consupmtion and labor(OLG loop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tolk=1e-3;              % Numerical tolerance for capital
tolc=1e-3;             % Numerical tolerance for energy consumption
tollab=1e-3;            % Numerical tolerance for labor
nq=10;                  % Max number of iterations

q=0;                    % Counter for iterations
K1=K0+10;
E1=E0+10;
N1=N0+10;

fprintf('\nComputing equilibrium price... \n');
while q<nq && (abs(K1-K0)>tolk ||abs(E1-E0)>tolc|| abs(N1-N0)>tollab)
   
    q=q+1;
    
    fprintf('\nIteration %g out of %g \n',q,nq);
    
     % Prices
    r0 = zeta*(((K0^zeta * N0^(1-zeta))^((phi-1)/phi)+(E0)^((phi-1)/phi))^(1/(phi-1)))*((K0^zeta)*(N0^(1-zeta)))^(-1/phi)-delta;
    w0 = (1-zeta)*(((((K0^zeta)*(N0^(1-zeta)))^((phi-1)/phi)+(E0)^((phi-1)/phi))^(1/phi-1)))*((K0^zeta)*(N0^(1-zeta)))^(-1/phi)*(K0/N0)^zeta;
   
    

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BACKWARD INDUCTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initialization
    v=zeros(nk,nk,J);                          % value function of agents(two dimensions:one is for capital number, one is for age)
    kapopt=ones(nk,J);                      % optimal savings of agents
    eneopt=ones(nk,J);                      % optimal energy consumption of agents
    % (store INDEX of k' in capital grid, not k' itself!)
    labopt=ones(nk,tW);                     % optimal labor supply   
    
        
    % Retired households
    
    % Last period utility
    cons=(1+r0)*kap-(Pe+tauC)*ene;   % last period consumption (vector!)
    ENE = transpose(ene);
    util=(((cons.^(gamma))*(ENE.^(1-gamma)))^(1-theta1))/(1-theta1);       % last period utility (vector!)(.:向量乘积）(l=0,since is the last period)
    v(:,J)=util;                            % last period indirect utility (vector!)
    
for j=J-1:-1:tW+1
    for ik0=1:nk
        vmin=neg;
        ik1=0;
 
        
    for ie0 = 1:nk
            vmin=neg;
            ie1=0;
            
            while ik1<nk && ie1<nk
            ik1=ik1+1;
            kap0=kap(ik0);
            kap1=kap(ik1);
            
           
            ie1=ie1+1;
            ene0=ene(ie0);
            ene1=ene(ie1);

     cons=(1+r0)*kap0-kap1-(Pe+tauC)*ene0;
      if cons<=0
                    util=neg;
                else
                    util=((cons^(gamma)* ene0^(1-gamma)))^(1-theta1)/(1-theta1);
      end
        
                  v0=util + beta*v(ik1,ie1,j+1);
          if v0>vmin
                    v(ik0,ie0,j)=v0;
                    kapopt(ik0,j)=ik1;
                    eneopt(ie0,1)=ie1;
                    vmin=v0;
          end
           end
            
     end
    end
end
  


        
    
    for j=tW:-1:1           % age
        for ik0=1:nk          % assets today
            
            % Initialize right-hand side of Bellman equation
            vmin=neg;
            ik1=0;
         for ie0 = 1:nk
            vmin=neg;
            ie1=0;
            
            % Loop over all k's in the capital grid to find the value,
            % which gives max of the right-hand side of Bellman equation
            
            while ik1 < nk && ie1<nk 	% assets tomorrow
                ik1=ik1+1;
                
                kap0=kap(ik0); % current asset holdings
                kap1=kap(ik1); % future asset holdings
                
                ie1=ie1+1;
                ene0=ene(ie0);
                ene1=ene(ie1);
                % Optimal labor supply
                % CHECK THIS!
                e(j)=e(tW-20,1);
                 lab=(gamma*(1-tauL)*e(j)*w0-(1-gamma)*((1+r0)*kap0-kap1))/((1-tauL)*e(j)*w0);
                
                % Check feasibility of labor supply
                if lab>1
                    lab=1;
                elseif lab<0
                    lab=0;
                end
            
                 
                
                % Instantaneous utility
                cons=(1+r0)*kap0+(1-tauL)*w0*e(j)*lab-kap1-(Pe+tauC)*ene0;%(e值要改）
                
                if cons<=0
                    util=neg;
                else
                    util=((cons^(gamma)* ene0^(1-gamma)))^(1-theta1)/(1-theta1)-1/(1+(1/theta2));
                end
            
        
                
                % Right-hand side of Bellman equation   
                v0=util +beta*v(ik1,ie1,j+1);
                
                % Store indirect utility, optimal saving and labor
                if v0>vmin
                    v(ik0,ie0,j)=v0;
                    kapopt(ik0,j)=ik1;
                    eneopt(ie0,j)=ie1;
                    labopt(ik0,ie0,j)=lab;
                    vmin=v0;
                end
            end
         end
         end
    end
    
            
        

           
   % Initializations
    kgen=zeros(J,1);        % capital supply k(j+1) for each generation j
    egen=zeros(J,1);
    labgen=zeros(tW,1);     % labor supply l(j) for each generation j

    % Use decision rules to iteratively fill in kgen and labgen
    ik0 = 1;                % starting capital of j = 1, kap(ik0) = 0(***)
    for j=1:J               % iterations over cohort
        % capital decision kp(k)(to find the optimal one)
        ik1 = kapopt(ik0, j);
        kgen(j) = kap(ik1);
        
        % energy decision e(k)
        ie1 = eneopt(ie0,j);
        egen(j) = ene(ie1);
        
        % labor decision l(k)
        if j<=tW
            labgen(j) = labopt(ik0, j);
        end

        % update k = kp(***)
        ik0 = ik1;
        ie0 =ie1;
    end

    K1 = kgen' * mass;              % dot product of vectors
    E1 = egen'* mass; 
    N1 = (labgen .* e)' * mass(1:tW);      % dot product of vectors

    % Update the guess on capital and labor    
    K0=0.9*K0+0.1*K1;
    E0=0.9*E0+0.1*E1;
    N0=0.9*N0+0.1*N1;

    % Display results
    disp('  capital     energy consumption     labor  ');
    disp([K0, E0, N0]);
    disp('deviation-capital   deviation-energy-consumption     deviation-labor       ');
    disp([abs(K1-K0),  abs(E1-E0), abs(N1-N0)]);
    
end


% Average hours worked
h = labgen'* mass(1:tW)/sum(mass(1:tW));

% Gini disposable income
income = zeros(J,1);
income(1:tW) = (1-tauL)*labgen .* e*w0 + r0*kgen(1:tW);
income(tW+1:end) = r0*kgen(tW+1:end);
gini_index = gini(mass,income);


% Display equilibrium results
% Output
Y=((K0^zeta*N0^(1-zeta))^((phi-1)/phi)+(E0)^((phi-1)/phi))^(phi/(phi-1));
% Income after labor taxation
incomeaftertax = (1-tauL)*w0;

disp('      K0         E0      N0       w         r     incomeaftertax    carbon taxation     Y     K/Y    h    Gini');
disp([K0, E0, N0, w0, r0,incomeaftertax, tauC, Y, K0/Y, h, gini_index]);




