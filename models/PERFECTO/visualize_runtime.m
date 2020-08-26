% PERFECTO: Prediction of Extended Response and Growth Functions for Estimating Chemotherapy Outcomes in Breast Cancer
% function to visualize network data at a given iteration in runtime
function id_maxv = visualize_runtime(populations, ~)
set(gcf, 'color', 'white');
% extract the max weight on each row (if multiple the first one)
id_maxv = zeros(populations(1).lsize, 1);
for idx = 1:populations(1).lsize
    [~, id_maxv(idx)] = max(populations(1).Wcross(idx, :));
end
HL = (rot90(populations(1).Wcross)); HL(HL<0)=0;
subplot(1, 2, 1);
hndl1 = imagesc(HL, [0, max(HL(:))]); box off; colorbar;
xlabel('neuron index'); ylabel('neuron index'); 
subplot(1, 2, 2);
hndl2 = surf(1:length(HL), 1:length(HL), HL);
colormap(jet)    % change color map
xlabel('neuron index'); ylabel('neuron index'); 
% refresh graphics
set(hndl1, 'CData', HL);
set(hndl2, 'CData', HL);
drawnow;
end