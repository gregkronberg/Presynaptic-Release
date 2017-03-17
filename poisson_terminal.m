function N = poisson_terminal(spikes,a,b,p0,ns,iter)

% simulate a single stochastic presynaptic terminal with a known spike train arriving at the
% synapse.  Vesicle docking and undocking are modeled as homogenous Poisson
% processes as in Zhang and Peskin 2015 PNAS
    % spikes = column vector of spike times
    % a = vesicle docking rate in (probability of vesicle docking per s)
    % b = vesicle undocking rate (probability of vesicle undocking per s)
    % p0 = probability of release for each vesicle when spike arrives
    % ns = number of vesicle sites per synapse
    % sims = number of iterations to run
    % N = number of vesicles released

N = zeros(iter,length(spikes));
parfor sim = 1:iter; % loop over simulations
    t = 0;
    pSite = a/(a + b); %steady state probability a site is occupied
    Ssite = rand(ns,1) < pSite; %state of sites (1 if occupied, 0 if empty)
    ap = 1; % spike index
    released = zeros(ns,length(spikes));
    while ap <= length(spikes); %loop through action potential times
        rates = (1-Ssite)*a + Ssite*b; % rate probability of state change at each site
        dt = -log(rand(ns,1))./rates; % poisson generated times to each event, scaled by rate probability
        [dtmin,i] = min(dt); % find the shortest time (i.e. the event most likely to happen first)
        if t+dtmin < spikes(ap); % if timestep remains with the current interval
            t = t + dtmin; % timestep to event occurance
            Ssite(i) = 1-Ssite(i); % change state of site
        else t = spikes(ap);% set time to spike arrival
            ri = rand(ns,1) < p0;% potential release sites
            released(:,ap) = ri.*Ssite;% if site is occupied, release
            Ssite = Ssite - released(:,ap); % change state at released sites
            ap = ap+1;% increase spike index
        end
    end
    N(sim,:) = sum(released,1); % store number of vesicles released
end
end