function boxplot_activity(PlotFolder,n_hour,dt_org,dt)
%boxplot_activity plots a boxplot figure for fish activities
%
% Inputs:
%      PlotFolder: folder of output figure
%      dt_org: original activity time inetrval in the resultfile, minutes
%      dt: required activity time interval to plot the activity boxplot
%
% Output:
%       Figure of the boxplot for each species
%

addpath('..\INCLUDE\export_fig')

% get the matlab result file
% via GUI 
[Fname, VidFolder] = uigetfile(fullfile('~' , ...
    '*.mat'), 'Select the matlab video file');
ResultFile = fullfile(VidFolder, Fname);

% load the matlab result file 
load(ResultFile,'PreferenceExperiment')

% Create string matrix for all fish names (IDs, such as COMP_XX or RYHN_XX)
Cases = fieldnames(PreferenceExperiment);
% Alternatevly choose some fishes
% Cases = {'COMP_31';'RHYN_01';'COMP_42';'RHYN_03'};

H = ['Hour_',num2str(n_hour,'%02.0f')];

% find how many fish in each specis
n_comp = 0;
n_rhyn = 0;
for i_case = 1:length(Cases)

    F = Cases{i_case};                                                     % fish name, COMP_XX or RHYN_XX
    
    if strcmp('COMP',F(1:4))
        n_comp = n_comp + 1;                                               % count up the COMP fish number
    elseif strcmp('RHYN',F(1:4))
        n_rhyn = n_rhyn + 1;                                               % count up the RHYN fish number
    else
        error('case is not COMP nor RHYN')
    end
    
end

% initialize variables
% time interval for activity
dt_r = dt_org/dt;                                                          % ratio of time interval
t_activity = 1/60:dt/60:n_hour;                                            % division by 60 is to convert minutes to hours
ActivSize = length(t_activity);
act1 = NaN(size(t_activity));
act2 = NaN(size(t_activity));
if n_comp > 0
    COMPActivBox = NaN(1,2*n_comp*ActivSize);
    comp_grp = COMPActivBox;
    comp_cases = strings(1,n_comp);
    comp_legndstr = strings(2*n_comp,1);
    i_comp = 0;
end
if n_rhyn > 0
    RHYNActivBox = NaN(1,2*n_rhyn*ActivSize);
    rhyn_grp = RHYNActivBox;
    rhyn_cases = strings(1,n_rhyn);
    rhyn_legndstr = strings(2*n_rhyn,1);
    i_rhyn = 0;
end

% loop for all fishes to get the required data
for i_case = 1:length(Cases)

    F = Cases{i_case};                                                     % name of current fish
    
    % 1.1 Rock
    temp1 = PreferenceExperiment.(F).RLSR.Analysis.(H).Activity;
    % 1.2 Sand:
    temp2 = PreferenceExperiment.(F).RRSL.Analysis.(H).Activity;
    for i=1:length(t_activity)
        ii = (i-1)*dt/dt_org+1;
        jj = ii+dt/dt_org-1;
        act1(i) = sum(temp1(ii:jj))*dt_r;
        act2(i) = sum(temp2(ii:jj))*dt_r;
    end
    
    % Data for COMP fishes
    if strcmp('COMP',F(1:4))
        i_comp = i_comp + 1; % to count up the COMP fish number
        % comp cases string
        comp_cases{i_comp} = F;
        % Box
        ii = 2*ActivSize*(i_comp-1)+1;
        jj = ii+ActivSize-1;
        ik = jj+1;
        jk = ik+ActivSize-1;
        COMPActivBox(1,ii:jj) = act1*100;
        comp_grp(1,ii:jj) = 2*(i_comp-1)+1;
        COMPActivBox(1,ik:jk) = act2*100;
        comp_grp(1,ik:jk) = 2*(i_comp);
        comp_legndstr(2*(i_comp-1)+1:2*i_comp,1) = num2str(i_comp);
        
        % Data for RHYN fishes
    else
        i_rhyn = i_rhyn + 1; % to count up the RHYN fish number
        % rhyn cases string
        rhyn_cases{i_rhyn} = F;
        % Box
        ii = 2*ActivSize*(i_rhyn-1)+1;
        jj = ii+ActivSize-1;
        ik = jj+1;
        jk = ik+ActivSize-1;
        RHYNActivBox(1,ii:jj) = act1*100;
        rhyn_grp(1,ii:jj) = 2*(i_rhyn-1)+1;
        RHYNActivBox(1,ik:jk) = act2*100;
        rhyn_grp(1,ik:jk) = 2*i_rhyn;
        rhyn_legndstr(2*(i_rhyn-1)+1:2*i_rhyn,1) = num2str(i_rhyn);

    end
end

% boxplot

% Comp
figure
set(gcf,'color','w');
boxplot(COMPActivBox,comp_grp,'Labels',convertStringsToChars(comp_legndstr))
title('Activity for C. compressirostris')
xlabel('C. compressirostris Indiveduals')
ylabel('Activity (%)')
ylim([0 100])

% RHYN
figure
set(gcf,'color','w');
boxplot(RHYNActivBox,rhyn_grp,'Labels',convertStringsToChars(rhyn_legndstr))
title('Activity for C. rhynchophorus')
xlabel('C. rhynchophorus Indiveduals')
ylabel('Activity (%)')
ylim([0 100])

end