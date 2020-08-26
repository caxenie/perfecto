% PERFECTO: Prediction of Extended Response and Growth Functions for Estimating Chemotherapy Outcomes in Breast Cancer
% crate the network composed of N_POP populations of
% N_NEURONS neurons implementing a SOM
% and init each struct weight and activity matrices
function populations = create_init_network(N_POP, N_NEURONS)
    wcross = rand(N_NEURONS, N_NEURONS);
    sigma_def = 0.045; % fixed to a half of the tuning curve normalized height
    for pop_idx = 1:N_POP
        populations(pop_idx) = struct(...
            'idx', pop_idx, ...
            'lsize', N_NEURONS, ...
            'Winput', zeros(N_NEURONS, 1),... % synaptogenesis
            's', sigma_def*ones(N_NEURONS, 1),...
            'Wcross', wcross./sum(wcross(:)), ...
            'a', zeros(N_NEURONS, 1));
    end
end