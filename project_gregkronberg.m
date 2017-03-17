% This program models a group of poisson presynaptic terminals that are
% stimulated by a focal extracellular electrode in a brain slice.
% A subset of these terminals is stochastically recruited by each 
% extracellular stimulus.  A train of stimuli is applied until a steady 
% state of vesicle release is achieved in the population.  A DC electric 
% field is then also applied, which modulates the membrane potential in axon 
% terminals.  This electric field is hypothesized to have two effects, which are compared
% in the model:  The first is to modulate the probability that a given
% stimulus will activate each axon, i.e. axon recruitment.  The second is 
% to modulate the probability of vesicle release, given that an action
% potential occured. Modulation of release probability appears to better 
% replicate experimental data.

%% Modulation of release probability with DC electric fields
clear all
close all
close all hidden

%% Synapse Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsyn = 50;% number of synapses with different parameters
ns = 10;% number of sites for each synapse
l = ones(ns,1);
k = 10; % sets scale for sigmoidal probability functions
stim = k*l*rand(1,nsyn); % stimulation strength felt at each synapse, proportional to proximity to stim electrode (arbitrary units)
pap = sigmf(stim,[k/4,k/2]);%probability that stimulus causes an action potential (sigmoidal function)
ca  = (k/3)*ones(ns,nsyn);% calcium entry per action potential (arbitrary units)  
p0 = sigmf(ca,[k/10,k/2]); % release probability per docked vesicle at time of action potential (sigmoidal function)
a = .1*l*ones(1,nsyn); % vesicle docking rate
b = .02*l*ones(1,nsyn); % vesicle undocking rate
    
%% Generate spike train
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
T = 15; % length of simulation in seconds
lambda = 10; % mean firing rate in Hz
spikes = (1:T*lambda)'/lambda;% column vector of spike times
        
%% DC electric field parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dc_on = .4*T;% onset of DC field
dc_off = .6*T;% offset of DC field
dstim_dc = .2*[1  0 -1   0]*k;% effective change in stimulation intensity with DC
dca_dc =   .2*[ 0 1   0 -1]*k;% effective change in calcium influx
color = {[0 0 1],[1 0 0],[0 0 1],[1 0 0]};
    
%% Run simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for q = 1:length(dstim_dc);
    pap_dc = sigmf(stim + dstim_dc(q),[k/4,k/2]);% change in action potential probability with DC
    p0_dc = sigmf(ca+dca_dc(q),[k/10,k/2]);%p0;%.25*(1-p0);% change in release probability with DC
    iter = 200;% number of interations
    N = zeros(iter,length(spikes),nsyn);
    parfor sim = 1:iter; % loop over simulations
        t = zeros(1,nsyn);% initialize time vector
        pSite = a./(a + b); %steady state probability a site is occupied
        Ssite = rand(ns,nsyn) < pSite; %state of sites (1 if occupied, 0 if empty)
        ap = 1; % spike index
        released = zeros(ns,nsyn,length(spikes));
        while ap <= length(spikes); %loop through action potential times
            rates = (1-Ssite(:)).*a(:) + Ssite(:).*b(:); % rate probability of state change at each site
            dt = -log(rand(ns*nsyn,1))./rates; % poisson generated times to each event, scaled by rate probability
            [dtmin,i] = min(dt); % find the shortest time for each synapse(i.e. the event most likely to happen first)
            if t+dtmin < spikes(ap); % if timestep remains with the current interval
                t = t + dtmin; % timestep to event occurance
                Ssite(i) = 1-Ssite(i); % change state of site
            else t = spikes(ap);% set time to spike arrival
                if t>=dc_on && t<=dc_off;% if dc field is on
                    ri = rand(ns,nsyn) < pap_dc.*p0_dc;% index sites with a release event
                else ri = rand(ns,nsyn) < pap.*p0;
                end
                released(:,:,ap) = ri.*Ssite;% if site is occupied, release
                Ssite = Ssite - released(:,:,ap); % change state at released sites
                ap = ap+1;% increase spike index
            end
        end
        released = permute(released,[1,3,2]);% rearrange release
        N(sim,:,:) = sum(released,1); % store number of vesicles released
    end
    
    %% Generate figures
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if rem(q,2)
        figure;hold on
    end
    base_i = ((dc_on/T)*lambda*T - 20):(dc_on/T)*lambda*T;
    N_sum = sum(N,3);
    N_norm = N_sum./(mean(N_sum(:,base_i),2)*ones(1,size(N_sum,2)));
    N_mean = mean(N_norm,1);
    N_sem = std(N_norm,0,1)/sqrt(size(N_norm,1));
    errorbar((spikes'),N_mean,N_sem,'.','Color',color{q});
    xSize = 15; ySize = 10;
    set(gcf,'PaperUnits','centimeters','Color', 'w')
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize],'Position',[0.5 0.5 xSize*50 ySize*50])
    set(gca,'TickDir','out','TickLength',[0.03 0.03],'box','off','FontSize',20); % make the tick sizes bigger and outside
    xlabel('time (s)');
    ylabel('Normalized vesicle release')
    xlim([4 12])
    if ~rem(q,2)
        legend('Recruitment','Release probability ')
    end
end
