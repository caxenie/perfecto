
## PERFECTO: Prediction of Extended Response and Growth Functions for Estimating Chemotherapy Outcomes in Breast Cancer

PERFECTO Codebase:

datasets - the experimental datasets (csv files) and their source, each in separate directories
models   - codebase to run and reproduce the experiments

Directory structure:

models/PERFECTO/.

create_init_network.m       - init PERFECTO network (SOM + HL)
error_std.m                 - error std calculation function
PERFECTO_core.m              - main script to run PERFECTO
model_rmse.m                - RMSE calculation function 
model_sse.m                 - SSE calculation function
parametrize_learning_law.m  - function to parametrize PERFECTO learning
present_tuning_curves.m     - function to visualize PERFECTO SOM tuning curves
randnum_gen.m               - weight initialization function
tumor_growth_model_fit.m    - function implementing ODE models
tumor_growth_models_eval.m  - main evaluation on PERFECTO runtime
visualize_results.m         - visualize PERFECTO output and internals
visualize_runtime.m         - visualize PERFECTO runtime


Usage: 

models/PERFECTO/perfecto_core.m - main function that runs PERFECTO and generates the runtime output file (mat file)
models/PERFECTO/tumor_growth_models_eval.m - evaluation and plotting function reading the PERFECTO runtime output file
