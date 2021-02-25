[~,msg] = unix('echo "$USER"'); % UNIX0platforme compliant only
if contains(msg,'omartin')
    path_oomao      = '/home/omartin/Projects/SIMULATIONS/OOMAO/oomao/'; % the oomao path
    path_workspace  = '/home/omartin/Projects/MAWSER/cnn_pyramid/'; % the simulation folder path
    path_save       = '/home/omartin/Projects/MAWSER/DICTS/';
elseif contains(msg,'matuser')
    path_oomao      = '/result/omartin/SIMULATIONS/OOMAO/mfiles/'; % the oomao path
    path_workspace  = '/result/omartin/Projects/CNN_WFS/cnn_pyramid/'; % the simulation folder path
    path_save       = '/result/omartin/Projects/CNN_WFS/DICTS/';
end

addpath(genpath(path_oomao),genpath(path_workspace));
