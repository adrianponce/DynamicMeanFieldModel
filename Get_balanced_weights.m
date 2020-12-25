%
%
% This prog. optimizes the strengh of the feedback inhibition of the FIC model 
% for varying global couplings (G)
% Saves the steady states and the feedback inhibition (J).
%
%
% For an isolated node, an input to the excitatory pool equal to Iie - be/ae = -0.026; 
% i.e., slightly inhibitory dominated, leads to a firing rate equal to 3.0631 Hz. 
% Hence, in the large-scale model of interconnected brain areas, 
% we aim to constraint in each brain area (i) the local feedback inhibitory weight Ji such 
% that Iie - be/aE = -0.026 is fulfilled (with a tolerance of ±0.005). 
% To achieve this, we apply following procedure: we simulate during 5000 steps 
% the system of stochastic differential DMF Equations and compute the averaged level of 
% the input to the local excitatory pool of each brain area,
% then we upregulate the corresponding local feedback inhibition Ji = Ji + delta;
% otherwise, we downregulate Ji = Ji - delta. 
% We recursively repeat this procedure until the constraint on the input
% to the local excitatory pool is fulfilled in all N brain areas.
%
% see:
% Deco et al. (2014) J Neurosci.
% http://www.jneurosci.org/content/34/23/7886.long
%
% Adrian Ponce-Alvarez
%--------------------------------------------------------------------------


clear all;

% Load connectome:
%--------------------------------
% load Human_66.mat C

load example_dti SC_asymmetric

N = size(SC_asymmetric,1);
ns = size(SC_asymmetric,3);

Cs = zeros(N,N,ns);
for j=1:ns
   C = SC_asymmetric(:,:,j); 
   C(1:N+1:end) = 0;
   Cs(:,:,j) = C/max( C(:) );
end

C=mean(Cs,3);

clear Cs


% Model's fixed parameters:
%----------------------------

dt=0.1;
tmax=10000;
tspan=0:dt:tmax;

taon=100;
taog=10;
gamma=0.641;
sigma=0.01;
JN=0.15;
J=ones(N,1);
I0=0.382; %%397;
Jexte=1.;
Jexti=0.7;
w=1.4;



% transfer functions:
% transfert function: excitatory
%--------------------------------------------------------------------------
ae=310;
be=125;
de=0.16;
He=@(x) (ae*x-be)./(1-exp(-de*(ae*x-be)));

% transfert function: inhibitory
%--------------------------------------------------------------------------
ai=615;
bi=177;
di=0.087;
Hi=@(x) (ai*x-bi)./(1-exp(-di*(ai*x-bi)));

% all tested global couplings:
wes  = .05:0.05:2; %.05:0.05:4.5; % warning: the range of wes depends on the conectome.
numW = length(wes);

% initialization:
%-------------------------

curr=zeros(tmax,N);
neuro_act=zeros(tmax,N);
delta=0.02*ones(N,1);

Se_init = zeros(N,numW);
Si_init = zeros(N,numW);

JI=zeros(N,numW);


kk=0;
for we=wes

    display(we)
    kk=kk+1;
    
    delta=0.02*ones(N,1);

    
%%% Balance (greedy algorithm)
% note that we used stochastic equations to estimate the JIs
% Doing that gives more stable solutions as the JIs for each node will be
% a function of the variance.

for k=1:5000
 sn=0.001*ones(N,1);
 sg=0.001*ones(N,1);
 nn=1;
 j=0;
 for i=2:1:length(tspan)
  xn=I0*Jexte+w*JN*sn+we*JN*C*sn-J.*sg;
  xg=I0*Jexti+JN*sn-sg;
  rn=feval(He,xn);
  rg=feval(Hi,xg);
  sn=sn+dt*(-sn/taon+(1-sn)*gamma.*rn/1000)+sqrt(dt)*sigma*randn(N,1);
  sn(sn>1) = 1; 
  sn(sn<0) = 0;             
  sg=sg+dt*(-sg/taog+rg/1000)+sqrt(dt)*sigma*randn(N,1);
  sg(sg>1) = 1;        
  sg(sg<0) = 0;
  j=j+1;
  if j==10
   curr(nn,:)=xn'-125/310;
   nn=nn+1;
   j=0;
  end
 end

currm=mean(curr,1);
 flag=0;
for n=1:1:N
  if abs(currm(n)+0.026)>0.005
   if currm(n)<-0.026 
    J(n)=J(n)-delta(n);
    delta(n)=delta(n)-0.001;
    if delta(n)<0.001;
       delta(n)=0.001;
    end
   else 
    J(n)=J(n)+delta(n);
   end
  else
   flag=flag+1;
  end
 end
 if flag==N
  break;
 end
end


 Se_init(:,kk)=sn; %Store steady states
 Si_init(:,kk)=sg; 
 JI(:,kk)=J; %Store feedback inhibition values
 
 % save results:
 
%  save Balanced_weights wes JI Se_init Si_init
 
end

save Benji_Balanced_weights wes JI Se_init Si_init



