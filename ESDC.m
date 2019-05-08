function [  ] = ESDC()
%opengl hardwarebasic

startup();

[input config] = input_processing();

disp('Starting analysis ...');
disp(' ');
disp('EP mass system scaling for reference points');
disp(' ');
data.reference = scale_EPsystem(input);
for i=1:size(input,2)
    disp(input(i).description);
    disp(data.reference(i));
end

% Simple analysis
data.analysis.F_const = analysis_F_const(input, data, config);
data.analysis.P_const = analysis_P_const(input, data, config);
%data.evolver_analysis = evolver_analysis(input, data, config)

% XML output
data = makestruct(input, data);
out.ESCD_output_data = data;

disp('Writing XML Output ')
struct2xml(out,'Output/ESDC_output')

disp('XML Output complete')
disp(' ')
disp('ESDC complete')
disp(' ')
disp(' ')
end

function output = makestruct(input, data)
disp('Preparing Output Data')
disp(' ')

output = struct();

for i=1:size(input,2)

    case_name = char(strcat('case',num2str(i)));

    output.(case_name).input = input(i);
    
    %structure calculation data for xml output
    data_fields = fieldnames(data);
    for k=1:size(data_fields ,1)
        if strcmp(data_fields(k),'reference')
           	for j=1:size(fieldnames(data.(char(data_fields(k)))),1)
                analysis_fields = fieldnames(data.(char(data_fields(k))));             
                output.(case_name).analysis.(char(data_fields(k))).(char(analysis_fields(j))) = data.(char(data_fields(k)))(i).(char(analysis_fields(j)));
            end
        else
            for j=1:size(fieldnames(data.(char(data_fields(k)))),1)
                analysis_fields = fieldnames(data.(char(data_fields(k))));             
                output.(case_name).analysis.(char(data_fields(k))).(char(analysis_fields(j))) = data.(char(data_fields(k))).(char(analysis_fields(j)))(i);
            end 
        end
    end


end

disp('Preparing Output Data Complete')
disp(' ')

end

function [analysis] = analysis_F_const(input, data, config)


%case specific walking of parameter space
disp('Simple analysis for c_e variation at constant thrust ...');
disp(' ');
%ToDo: Add plot for P_electric
%ToDo: Add plot for burn time


analysis = input ; 

for i=1:size(input,2)  % for each concept
    
   analysis(i).description = char([analysis(i).description ' and F = '  num2str(analysis(i).F) ' N']);
        c_e = linspace(config.(analysis(i).thruster_type).c_e_min, config.(analysis(i).thruster_type).c_e_max, 1000);
         % Todo: adapt resolution = 1000 through config here
        % Calculate data from c_e span    
        analysis(i).c_e = c_e;
        analysis(i).m_dot= analysis(i).F./c_e;
        analysis(i).P_jet = 1/2.*analysis(i).m_dot.*c_e.^2;
        
        analysis(i).PPU_P = analysis(i).P_jet./( analysis(i).PPU_eff*analysis(i).eff);

        analysis(i).result = scale_EPsystem(analysis(i));
        disp(' ')
        disp('Min achieveable \mu_{EP}');
        disp(analysis(i).result.mu_EP_min);
        disp('at');
        disp(strcat('c_e = ',num2str(analysis(i).result.mu_EP_min_c_e),' m/s'));
        disp(' ');

    %Output
    % Plot Component Sysem Mass fraction of components normed to total sys mass
        if config.plots_componentwise == 1
                makePlotNormalizedSystemFunction(analysis(i));
        end
end

% Plot comparison of system mass fractions of all input concepts
    if config.plot_systemcompare == 1  
            makePlotFullSystemComparison(input, data, analysis,config, ' constant thrust');        
    end
    
    
    
    makePlotFullSystemPowerComparison(analysis, ' constant thrust'); 
    makePlotFullSystemBurnTimeComparison(analysis,' constant thrust');

end

function [analysis] = analysis_P_const(input, data, config)
disp('Simple analysis for c_e variation at constant power ...');
disp(' ');
%ToDo: Add plot for resulting F
%ToDo: Add plot for burn time

analysis = input ; 




for i=1:size(input,2)  % for each concept
    analysis(i).description= char([analysis(i).description ' and P = ' num2str(analysis(i).P_propulsion) ' W']);

        c_e = linspace(config.(analysis(i).thruster_type).c_e_min, config.(analysis(i).thruster_type).c_e_max, 1000);
         % Todo: adapt resolution = 1000 through config here

        % Calculate data from c_e span    
      analysis(i).c_e = c_e;
      analysis(i).P_jet = analysis(i).P_propulsion/(analysis(i).eff*analysis(i).PPU_eff).*ones(1,size(c_e,2));


   %   analysis(i).m_dot= 2.*analysis(i).P_jet./(c_e.^2);
      analysis(i).PPU_P = analysis(i).P_propulsion.*ones(1,size(c_e,2));
      
      analysis(i).result = scale_EPsystem(analysis(i));

      

      disp('Min achieveable \mu_{EP}');
      disp(analysis(i).result.mu_EP_min);
      disp('at');
      disp(strcat('c_e = ',num2str(analysis(i).result.mu_EP_min_c_e),' m/s'));
      disp(' ');
      
%ToDo: Constant electrical power supply


        if config.plots_componentwise == 1
                makePlotNormalizedSystemFunction(analysis(i));
        end
        
        %add plot thrust progression 

end
% Plot comparison of system mass fractions of all input concepts
    if config.plot_systemcompare == 1  
            makePlotFullSystemComparison(input, data, analysis,config, ' constant power supply');        
    end
    
    makePlotFullSystemThrustComparison(analysis,' constant power supply'); 
    makePlotFullSystemBurnTimeComparison(analysis,' constant power supply');

end

%Evolver 

function [data] = evolver_analysis(input, data, config) % todo more configs with respect to evolver    
disp('Evolver analysis starting...');
disp(' ');

landscape =  make3DFitnesslandscape(input, data, config);


%make population
population = struct();
[population.initial] = makeInitialPopulation(input, config);


%make fitness landscapes todo here
%TODO landscape each case and resulting min landscape


%evolve
evolutionaryConverger(population, input, data, config);

end

function [initial_pop] = makeInitialPopulation(input, config)

initial_pop = struct();
n_cases = size(input,2);

    for i=1:config.n_seeds
        fieldname = strcat('seed',num2str(i));
        initial_pop.(fieldname).case = randi(n_cases);

        c_e_min =   config.(char(input(initial_pop.(fieldname).case).thruster_type)).c_e_min;
        c_e_max =   config.(char(input(initial_pop.(fieldname).case).thruster_type)).c_e_max;
        F_min   =   config.(char(input(initial_pop.(fieldname).case).thruster_type)).F_min;
        F_max   =   config.(char(input(initial_pop.(fieldname).case).thruster_type)).F_max;

        initial_pop.(fieldname).c_e  = c_e_min + (c_e_min+c_e_max)*rand;
        initial_pop.(fieldname).F    = F_min + (F_min+F_max)*rand;
    end
    
end

function [population_history] = evolutionaryConverger(inital_population, m_0, eta_PPU, m_struct, I_tot, x_mesh, y_mesh, z_mesh, colors)

%t_start = now;
%datetime('now');
convergence = 0;
% todo self convergence later
%eps_total_converge = 0.001;
population_history = {};

%video stuff here
%video=VideoWriter('Evolver_testvid')
%open(video)
%close(video)

current_population = ones(size(inital_population,1),4);
new_population = ones(size(inital_population,1),4);
current_mu = ones(size(inital_population,1),1);
    
 for i=1:size(current_population,1)
          F         = inital_population(i, 1);
          c_e       = inital_population(i, 2);

          prop_type = inital_population(i, 3);
          if prop_type==1
              eta_PPU_current = eta_PPU(1);
              I_tot_in = I_tot(1);
              m_struct_in = m_struct(1);
              
          else
              eta_PPU_current = eta_PPU(2);        
              I_tot_in = I_tot(2);
              m_struct_in = m_struct(2);
          end
          current_mu(i) =   calculateMassFraction(F, c_e, m_0, eta_PPU_current, prop_type,m_struct_in, I_tot_in,1);
 end
 
 current_population = [inital_population current_mu];
 population_history{end+1}= current_population;


generation = 0;

while convergence == 0
    generation = generation +1;
%     gen_frame = figure('Name', 'Converger Full History');
%     mesh(x_mesh,y_mesh, z_mesh, colors, 'FaceAlpha', 0.1, 'EdgeAlpha', 0.5);
%     view(-142,30)
%     xlim([0.01 0.1])
%     ylim([0 max(get_c_e_span())])
%     zlim([0 1])
%     xlabel('F / N');
%     ylabel('c_e / m/s');
%     zlabel('\mu_{EP} / - ');
%     hold on
    for i=1:size(current_population,1)
        %mutator
        F           = mutate_F(current_population(i, 1));
        c_e         = mutate_c_e(current_population(i, 2));
        prop_type   = mutate_prop(current_population(i, 3));
           if prop_type==1
              eta_PPU_current = eta_PPU(1);
              I_tot_in = I_tot(1);
              m_struct_in = m_struct(1);
              
          else
              eta_PPU_current = eta_PPU(2);        
              I_tot_in = I_tot(2);
              m_struct_in = m_struct(2);
          end
          new_mu=   calculateMassFraction(F, c_e, m_0, eta_PPU_current, prop_type,m_struct_in, I_tot_in,1);

   % plot3([current_population(i, 1)], [current_population(i, 2)] ,[new_mu],'rx','LineWidth',2);
     %  set(gcf,'PaperPositionMode','auto')

               
        %selector
        if new_mu <= current_mu(i)
            new_population(i,:) = [F c_e prop_type new_mu];
            current_mu(i) = new_mu;
            
        else
            new_population(i,:) = current_population(i,:);
        end
    end
%     hold off
%     framename= sprintf('frame%03d', generation);
%     print(gen_frame, framename,'-dpng')
%     close all
    
population_history{end+1}=new_population;
    
current_population = new_population;
%have convergence when epsilon final is reached.

%iterate convergence distance once hyperspace velocity is low
%     if generation == 10
%     title = strcat('Mass fraction EP System with ', num2str(size(inital_population,1)),' inital seed points for',num2str(generation),' generations');
%     figure('Name', title);
%     mesh(x_mesh,y_mesh, z_mesh, colors, 'FaceAlpha', 0.1, 'EdgeAlpha', 0.5);
%     view(-142,30)
%     hold on
% 
% for i=1:size(population_history{1},1)
% x=[];
% y=[];
% z=[];
%     for k=1:numel(population_history)
%       pop=  population_history{k};
% 
%         x = [x pop(i,1)];
%         y = [y pop(i,2)];
%         z = [z pop(i,4)];   
%     end
%     plot3(x,y,z,'color',rand(1,3));
% end
%     xlabel('F / N');
%     ylabel('c_e / m/s');
%     zlabel('\mu_{EP} / - ');
%     
%     hold off
%     end
    
%     if generation == 100
%    title = strcat('Mass fraction EP System with ', num2str(size(inital_population,1)),' inital seed points for',num2str(generation),' generations');
%     figure('Name', title);
%     mesh(x_mesh,y_mesh, z_mesh, colors, 'FaceAlpha', 0.1, 'EdgeAlpha', 0.5);
%     view(-142,30)
%     hold on
% 
% 
% 
%     for i=1:size(population_history{1},1)
%     x=[];
%     y=[];
%     z=[];
%         for k=1:numel(population_history)
%           pop=  population_history{k};
% 
%             x = [x pop(i,1)];
%             y = [y pop(i,2)];
%             z = [z pop(i,4)];   
%         end
%     plot3(x,y,z,'color',rand(1,3));
%     end
%     xlabel('F / N');
%     ylabel('c_e / m/s');
%     zlabel('\mu_{EP} / - ');
%     
%     hold off     
%     end
    
    if generation == 1000
%         hold off
%         title = strcat('Mass fraction EP System with ', num2str(size(inital_population,1)),' inital seed points for',num2str(generation),' generations');
%     figure('Name', title);
%     mesh(x_mesh,y_mesh, z_mesh, colors, 'FaceAlpha', 0.1, 'EdgeAlpha', 0.5);
%     view(-142,30)
%     hold on



% for i=1:size(population_history{1},1)
% x=[];
% y=[];
% z=[];
%     for k=1:numel(population_history)
%       pop=  population_history{k};
% 
%         x = [x pop(i,1)];
%         y = [y pop(i,2)];
%         z = [z pop(i,4)];   
%     end
%     plot3(x,y,z,'color',rand(1,3));
% 
% 
% end
%     xlabel('F / N');
%     ylabel('c_e / m/s');
%     zlabel('\mu_{EP} / - ');
%     
%     hold off
    convergence =1;
    end
end

% [global_min  glob_min_index]= min(z);
% mu_min = global_min;
% F_min = x(glob_min_index);
% c_e_min = y(glob_min_index);

%close(video)
t_end = now;
datetime('now');
dt = t_start - t_end;
end

function [c_e_new] = mutate_c_e(c_e_old)
c_e_span = get_c_e_span();

resolution_const=0.0001; % here potential for improvement
diff = (max(c_e_span)-min(c_e_span));
direction = (-1)^randi(2);
increment = randi(100); % potential for improvement

c_e_new = c_e_old+direction*increment*diff*resolution_const;

%mutation definition here

if c_e_new> max(c_e_span)
    c_e_new = max(c_e_span);
end
if c_e_new< min(c_e_span)
    c_e_new = min(c_e_span);
end
end

function [prop_new] = mutate_prop(prop_old)
%only mutate between applied propellants - no increments
 prop_new = randi(max(get_prop_span()));
end

function [F_new] = mutate_F(F_old)
F_span = get_F_span();
resolution_const=0.001; % here potential for improvement
diff = (max(F_span)-min(F_span));
direction = (-1)^randi(2);
increment = randi(10); %10 potential for improvement

F_new = F_old+direction*increment*diff*resolution_const;

%mutation definition here

if F_new> max(F_span)
    F_new = max(F_span);
end
if F_new< min(F_span)
    F_new = min(F_span);
end

end

%Scaling 

function [data] = scale_EPsystem(input)


for i = 1:size(input,2)
    
    m = linspace(0,0,size(input(1).c_e,2));    
    
    data(i).PPU_mass = m_scale_PPU(input(i).PPU_P,input(i).PPU_eff, input(i).thruster_type);
    m = m +data(i).PPU_mass;

    data(i).PV_mass = m_scale_SolarPanel(input(i).PPU_P);
    m = m + data(i).PV_mass;

    data(i).propellant_mass = m_scale_propellant(input(i).sat_mass, input(i).delta_v, input(i).c_e);
    m = m + data(i).propellant_mass;

    data(i).tank_mass = m_scale_tank(data(i).propellant_mass, input(i).propellant);
    m = m + data(i).tank_mass;

    data(i).thruster_mass = m_scale_thruster(input(i).PPU_P, input(i).thruster_type); 
    m = m + data(i).thruster_mass;
   
    
    data(i).structure_mass = m_scale_structure(input(i).struct_mass);
    m = m + data(i).structure_mass;

    EP_sys_mass= m;
        for j=1:size(EP_sys_mass,2)
            if EP_sys_mass(j)>input(i).sat_mass
                EP_sys_mass(j) = input(i).sat_mass+1;
            end
        end
    data(i).EP_sys_mass = EP_sys_mass;

    mu_EP = data(i).EP_sys_mass./input(i).sat_mass;
        for j=1:size(mu_EP,2)
            if mu_EP(j)>1
                mu_EP(j) = 1.01;
            end
        end
    data(i).mu_EP =mu_EP;
    
    if size(input(i).c_e,2)>1
        [data(i).mu_EP_min, mu_min_index] = min(data(i).mu_EP);
        data(i).mu_EP_min_c_e = input(i).c_e(mu_min_index);
        PPU = data(i).PPU_mass(mu_min_index);
        PV  =data(i).PV_mass(mu_min_index);
        tank = data(i).tank_mass(mu_min_index);
        thruster = data(i).thruster_mass(mu_min_index);
        prop = data(i).propellant_mass(mu_min_index);
    end

    
end


end

function [m_PPU_out] = m_scale_PPU(P_in,eta_PPU,type)

    if strcmp(type,'Arcjet')
% Ref:BB28 report for scaling of PPU mass
    m_0     = 1.011;
    slope   = 2.465/1000;
    
    elseif strcmp(type,'GridIonThruster')
            m_0     = 1.5;       % DCIU % AIAA-2006-5 162  Mission Benefits of Gridded Ion and Hall Thruster Hybrid Propulsion Systems
            slope   = 2.5/1000;
        
    else
        disp('Unknown thruster - no PPU scaling available in m_scale_PPU')
    end
    
    m_PPU_out = m_0+slope.*P_in./eta_PPU;
end

function [m_Panel_out] = m_scale_SolarPanel(P_out)
%currently scaled for requested output power
    m_0     = 0;
    slope   = 1/64.242; % W/kg Ref:http://dhvtechnology.com/wp-content/uploads/2017/07/Datasheet-Julio-v1-front-back.pdf

    m_Panel_out = m_0+ slope.*P_out;

end

function [m_Tank_out] = m_scale_tank(m_prop, type)


%https://www.orbitalatk.com/commerce/Data_Sheet_Index_Diaphragm-VOL.aspx

m_0_set =   [];
slope_set = [];
rho_p_set = [];
%https://www.orbitalatk.com/commerce/Data_Sheet_Index_Diaphragm-VOL.aspx
for i=1:size(m_prop,2)

    if strcmp(type,'He') % Helium
        m_0     = 2.77;
        slope   = 250.8;
        rho_p   = 43.14;    % at 300 bar 20°C https://www.wolframalpha.com/input/?i=density+of+helium+at+300+bar
    end

    if strcmp(type,'NH3') % Ammonia
        m_0     = 1.45;
        slope   = 282.; 
        rho_p   = 681.9;
    end
    
    if strcmp(type,'Xe') % Xenon at 20°C https://www.wolframalpha.com/input/?i=xenon+density+at+300+bar+at+20+%C2%B0C
        m_0     = 2.77;
        slope   = 250.8;
        rho_p   = 2300;
    end
    m_0_set     =   [m_0_set m_0];
    slope_set   =   [slope_set slope];
    rho_p_set   =   [rho_p_set rho_p];

end
    m_Tank_out = m_0_set+ slope_set.*(m_prop./rho_p_set).^(3/2);

end

function [m_Thruster_out] = m_scale_thruster(P_in, type)

m_0_set = [];
slope_set = [];

for i=1:size(P_in,2)
    
    if strcmp(type,'Arcjet')
    
        if P_in(i) <= 300
            m_0     =0.3;        % Velarc mass
            slope   =0;        % constant for P <300 W
        end
        if ((P_in(i)>= 300) && (P_in(i) < 1500))
            m_0     = 0.2;        % ATOS/ARTUS 0.7 kg at 1500 W
            slope   = 0.00033333333333333333;

        end
        if ((P_in(i) >= 1500) && (P_in(i) < 10000))
            m_0     = 0.4882352941176471;        % MARC
            slope   = 0.00014117647058823529;
        end

        if P_in(i) >= 10000
            m_0     = 0.0196581196581199205;        % HIPARC
            slope   = 0.00018803418803418803;  
        end
    m_0_set = [m_0_set m_0];
    slope_set = [slope_set slope];
    elseif strcmp(type,'GridIonThruster')       % http://www.space-propulsion.com/brochures/electric-propulsion/electric-propulsion-thrusters.pdf
            
        if P_in(i) <= 50                            % RIT μX
            m_0     =0.44;    
            slope   =0;             
        end
        
        if ((P_in(i)>= 50) && (P_in(i) < 2000))     %RIT 10 EVO
            m_0     =0.027;    
            slope   =0.0034;       
        end
        
        if P_in(i) >= 2000                          % RIT 2X
            m_0     =0.027;    
            slope   =0.0034;    
        end
    m_0_set = [m_0_set m_0];
    slope_set = [slope_set slope];
    else
        disp('This type of thruster not yet implemented');
    end
end
    m_Thruster_out = m_0_set + slope_set.*P_in;
end

function [m_propellant_out] = m_scale_propellant(mass, dv, c_e)
    I_tot = c_e.*mass.*(1-exp(-dv./c_e));
    m_propellant_out = I_tot./c_e;
end

function [m_struct_out] = m_scale_structure(m)
% elaborate further?

m_struct_out= m;
end

%Input

function [] = startup()
clc;
close all;
disp('Start ESDC solver tool')
disp(' ')
disp('Institute of Space Systems')
disp('Author: Manfred Ehresmann')
disp(' ')
disp('Loading Files ...')
disp(' ')
end

function [input_struct] = read_input_mission_parameter()
 disp('Reading Mission Paramater Input File');
% read m_0

input_struct = struct();

mission_input = xmlread('Input/ESDC_Input.xml');
satellite_parameters = mission_input.getDocumentElement;


%cases
allcaseitems = satellite_parameters.getElementsByTagName('case');


%walk xml
    for k=0:allcaseitems.getLength-1
        thiscase= allcaseitems.item(k);
        %total_impulse
        Item = thiscase.getElementsByTagName('totalmass');
        input_struct(k+1).sat_mass       = str2double(Item.item(0).getFirstChild.getData);
        %total_impulse
        Item = thiscase.getElementsByTagName('deltav');
        input_struct(k+1).delta_v     = str2double(Item.item(0).getFirstChild.getData);
        %thruster
        Item = thiscase.getElementsByTagName('thruster');
        input_struct(k+1).thruster_type  = char(Item.item(0).getFirstChild.getData);
        %propellamnt
        Item = thiscase.getElementsByTagName('propellant');
        input_struct(k+1).propellant = char(Item.item(0).getFirstChild.getData); 
        Item = thiscase.getElementsByTagName('P_propulsion');
        input_struct(k+1).P_propulsion = str2double(Item.item(0).getFirstChild.getData); 
    end

 disp('Success');
 disp(' ');
end

function [configuration] = read_input_simulation_parameter()
 disp('Reading Simulation Parameter Input File')
 
 configuration = struct();
 
 sim_input = xmlread('Input/ESDC_Simulation_parameters.xml');
 sim_parameters = sim_input.getDocumentElement;
 
 prop_sys = {'Arcjet' ; 'GridIonThruster'}; %todo automate generation of this list
 

 for i=1:size(prop_sys,1)
    name = char(prop_sys(i));
    parameter_type  = sim_parameters.getElementsByTagName(name);
    parameter       = parameter_type.item(0).getElementsByTagName('c_e');
    parameter_min   = parameter.item(0).getElementsByTagName('min');
    configuration.(name).c_e_min         = str2double(parameter_min.item(0).getFirstChild.getData);
    parameter_max   = parameter.item(0).getElementsByTagName('max');
    configuration.(name).c_e_max         = str2double(parameter_max.item(0).getFirstChild.getData);
    parameter       = parameter_type.item(0).getElementsByTagName('F');
    parameter_min   = parameter.item(0).getElementsByTagName('min');
    configuration.(name).F_min           = str2double(parameter_min.item(0).getFirstChild.getData);
    parameter_max   = parameter.item(0).getElementsByTagName('max');
    configuration.(name).F_max           = str2double(parameter_max.item(0).getFirstChild.getData);

 end
     
    %System plot
    sysplot= sim_parameters.getElementsByTagName('system');
    configuration.plot_systemcompare= str2double(sysplot.item(0).getFirstChild.getData);
    
    %Component plot
    componentplot= sim_parameters.getElementsByTagName('components');
    configuration.plots_componentwise= str2double(componentplot.item(0).getFirstChild.getData);
    
    %margin 
    margin_data = sim_parameters.getElementsByTagName('margin');
    configuration.margin = str2double(margin_data.item(0).getFirstChild.getData);
 
    %evolver
    seed_points = sim_parameters.getElementsByTagName('seed_points');
    configuration.n_seeds = str2double(seed_points.item(0).getFirstChild.getData);
    
    landscape_res = sim_parameters.getElementsByTagName('fitness_landscape_resolution');
    configuration.res_landscape = str2double(landscape_res.item(0).getFirstChild.getData);
 
 disp('Success')
 disp(' ')
end

function [reference_data]  = read_reference_data()
    disp('Reading Reference Data Input File');
    
    reference_data = struct();
    
    ref_data_input = xmlread('Database/ESDC_Reference_Data.xml');
    ref_data = ref_data_input.getDocumentElement;
    
    
    allcaseitems = ref_data.getElementsByTagName('case');
    

    for k=0:allcaseitems.getLength-1
        thiscase = allcaseitems.item(k);
       
        type_xml = thiscase.getElementsByTagName('type');
        reference_data(k+1).thruster_type = char(type_xml.item(0).getFirstChild.getData);
        
        name_xml = thiscase.getElementsByTagName('name');
        reference_data(k+1).name = char(name_xml.item(0).getFirstChild.getData);
        
        prop_xml = thiscase.getElementsByTagName('propellant');
        reference_data(k+1).propellant = char(prop_xml.item(0).getFirstChild.getData);
        
        thrust_xml = thiscase.getElementsByTagName('thrust');
        reference_data(k+1).F = str2double(thrust_xml.item(0).getFirstChild.getData);
        
        c_e_xml = thiscase.getElementsByTagName('c_e');
        reference_data(k+1).c_e = str2double(c_e_xml.item(0).getFirstChild.getData); 
        
        m_dot_xml = thiscase.getElementsByTagName('massflow');
        reference_data(k+1).m_dot = str2double(m_dot_xml.item(0).getFirstChild.getData); 
        
        mass_xml = thiscase.getElementsByTagName('mass');
        reference_data(k+1).mass = str2double(mass_xml.item(0).getFirstChild.getData); 
        
        eff_thruster_xml = thiscase.getElementsByTagName('efficiency');
        reference_data(k+1).eff = str2double(eff_thruster_xml.item(0).getFirstChild.getData); 
         
        m_struct_xml = thiscase.getElementsByTagName('structure_mass');
        reference_data(k+1).mass_struct = str2double(m_struct_xml.item(0).getFirstChild.getData); 
        
        PPU_P_xml = thiscase.getElementsByTagName('PPU_P');
        reference_data(k+1).PPU_P = str2double(PPU_P_xml.item(0).getFirstChild.getData);  
        
        PPU_eff_xml = thiscase.getElementsByTagName('PPU_efficiency');
        reference_data(k+1).PPU_eff = str2double(PPU_eff_xml.item(0).getFirstChild.getData);       
    end
        

 disp('Success');
 disp(' ');
end 

function [input configuration]= input_processing()
 %InputReading
    
    % Mission Parameters
    [input] =  read_input_mission_parameter();
   
    % Database     
    [ref_data] =  read_reference_data();
    
    %[thruster_type_ref, name, propellant_ref, F, c_e, m_dot, mass_ref, eff, mass_struct_ref, PPU_P, PPU_eff] 
    
    for i=1:size(input,2)

        for j=1:size(ref_data,2)
            
            
            if and(strcmp(input(i).thruster_type,ref_data(j).thruster_type),strcmp(input(i).propellant,ref_data(j).propellant))
                case_index=j;
            end
        end
        
    input(i).description=char(strcat({'Case: '},input(i).thruster_type, {' with '} , input(i).propellant, {[char(10) ' for ']}, num2str(input(i).delta_v), {' m per s'})); 
                      input                                                                                         
    input(i).name       = ref_data(case_index).name;
    input(i).F          = ref_data(case_index).F;
    input(i).c_e        = ref_data(case_index).c_e;
    input(i).m_dot      = input(i).F/input(i).c_e;
    input(i).eff        = ref_data(case_index).eff;
    input(i).thruster_mass = ref_data(case_index).mass;
    input(i).struct_mass= ref_data(case_index).mass_struct;
    input(i).PPU_P      = ref_data(case_index).PPU_P;
    input(i).PPU_eff    = ref_data(case_index).PPU_eff;
    input(i).P_jet      = 1/2*input(i).m_dot *input(i).c_e^2;
    input(i).propellant_mass = m_scale_propellant(input(i).sat_mass , input(i).delta_v ,  input(i).c_e);

    end
    
        % Simulation Parameters
    configuration = read_input_simulation_parameter();
    if configuration.margin >=0
             configuration.ConsiderMargins = 1;
         else
             configuration.ConsiderMargins = 0;
    end
    
    disp('Simulation configuration:');
    disp(configuration);
    disp(' ');
    
    disp('Input cases:');
    for i=1:size(input,2)
        disp(input(i));
        disp(' ');
    end
    
    disp(' ');
    disp('Input Reading complete');
    disp(' ');
end

%Output plots

function [] = makePlotNormalizedSystemFunction(input)
disp(' Output production: Normalized System Component Graphs');
disp(' ')

    fig_handle = figure('Name',strcat('Component comparison ',input.description));
    hold on
        plot(input.c_e,input.result.EP_sys_mass./input.result.EP_sys_mass,'-k','LineWidth',1.5);
        plot(input.c_e,input.result.PPU_mass./input.result.EP_sys_mass,'-b','LineWidth',1.5);
        plot(input.c_e,input.result.PV_mass./input.result.EP_sys_mass,'-r','LineWidth',1.5);
        plot(input.c_e,input.result.tank_mass./input.result.EP_sys_mass,'-c','LineWidth',1.5);
        plot(input.c_e,input.result.thruster_mass./input.result.EP_sys_mass, 'color',[0 .5 .5],'LineWidth',1.5);
        plot(input.c_e,input.result.structure_mass./input.result.EP_sys_mass,'-y','LineWidth',1.5);
        plot(input.c_e,input.result.propellant_mass./input.result.EP_sys_mass,'-m','LineWidth',1.5);
        plot([input.result.mu_EP_min_c_e input.result.mu_EP_min_c_e],[0 1],'-dk');
    xlabel('c_e / m/s')
    ylabel('\mu_{Component} / -')
    ylim([0 1])
    %todo plot2 file 
    
    legend('System','PPU','Solar Panel','Tank','Thruster','Structure','Propellant','Minimum');

    saveas(fig_handle,strcat('Output/','ESDC Component Comparison Plot ',strrep(input.description,char(10),''),'.png'),'png');

    disp(strcat('Output complete: Componentwise comparison',input.description));
        hold off
end

function [] = makePlotFullSystemComparison(input, data, analysis,config, case_type)
disp(' ')
disp(' Output production: Full System Comparison Graph');
disp(' ')
plotLegend = [];

    if config.ConsiderMargins == 1
        analysis = plot_margins(analysis, config);
    end

    fig_handle = figure('Name',strcat('Full Case Comparison ', case_type));
    hold on;
    for i=1:size(analysis,2)
        plot(analysis(i).c_e,analysis(i).result.mu_EP,'color',rand(1,3),'LineWidth',1.5 );
        plotLegend=[plotLegend cellstr(analysis(i).description)];
    end
    
    plotLegend =[ plotLegend 'Minima'];
    plotLegend =[ plotLegend 'Current Design'];
    
        if config.ConsiderMargins == 1
            plotLegend=[ plotLegend strcat(num2str(config.margin), ' % Margin')];
        end

    for i=1:size(analysis,2)
        plot(analysis(i).result.mu_EP_min_c_e,analysis(i).result.mu_EP_min,'dr');
        plot(input(i).c_e,data.reference(i).mu_EP,'xr');
            if config.ConsiderMargins == 1
               plot(analysis(i).mu_min_margin_c_e, analysis(i).mu_min_margin, '--k','LineWidth',1.5);
            end
    end
    
    legend(plotLegend);
    xlabel('c_e / m/s')
    ylabel('\mu_{EPropSys} / -')
    ylim([0 1])
    
    saveas(fig_handle,strcat('Output/','ESDC Full Concept Comparison ', case_type,'.png'),'png')
    disp(strcat('Output complete: Full Concept Comparison Plot'));
end

function [] = makePlotFullSystemBurnTimeComparison(analysis,case_type)
disp(' ')
disp(' Output production: Full System Burn Time Comparison Graph');
disp(' ')
plotLegend = [];

    fig_handle = figure('Name',strcat('Full Case  Burn Time Comparison ', case_type));
    hold on;
    for i=1:size(analysis,2)
        t_burn = analysis(i).result.propellant_mass./analysis(i).m_dot;
        t_burn = t_burn./(60*60*24);
        plot(analysis(i).c_e,t_burn,'color',rand(1,3),'LineWidth',1.5 );
        plotLegend=[plotLegend cellstr(analysis(i).description)];
    end
    hold off

    legend(plotLegend);
    xlabel('c_e / m/s')
    ylabel('t_{burn} / days')
    
    saveas(fig_handle,strcat('Output/','ESDC Full Concept Burn Time Comparison ', case_type,'.png'),'png')
    disp(strcat('Output complete: Full Concept Burn Time Comparison Plot', case_type));
end

function [] = makePlotFullSystemPowerComparison(analysis,case_type)
disp(' ')
disp(' Output production: Full System Power Comparison Graph');
disp(' ')
plotLegend = [];

    fig_handle = figure('Name',strcat('Full Case Power Comparison ', case_type));
    hold on;
    for i=1:size(analysis,2)
        plot(analysis(i).c_e,analysis(i).P_jet,'color',rand(1,3),'LineWidth',1.5 );
        plotLegend=[plotLegend cellstr(analysis(i).description)];
    end
    hold off

    legend(plotLegend);
    xlabel('c_e / m/s')
    ylabel('P_{jet} / W')
    
    saveas(fig_handle,strcat('Output/','ESDC Full Concept Power Comparison ', case_type,'.png'),'png')
    disp(strcat('Output complete: Full Concept Power Comparison Plot'));
end

function [] = makePlotFullSystemThrustComparison(analysis,case_type)
disp(' ')
disp(' Output production: Full System Comparison Graph');
disp(' ')
plotLegend = [];

    fig_handle = figure('Name',strcat('Full Case Comparison ', case_type));
    hold on;
    for i=1:size(analysis,2)
        F = 2*analysis(i).P_jet./analysis(i).c_e;
        plot(analysis(i).c_e,F,'color',rand(1,3),'LineWidth',1.5 );
        plotLegend=[plotLegend cellstr(analysis(i).description)];
    end
    hold off

    legend(plotLegend);
    xlabel('c_e / m/s')
    ylabel('F / N')

    saveas(fig_handle,strcat('Output/','ESDC Full Concept Comparison ', case_type,'.png'),'png')
    disp(strcat('Output complete: Full Concept Comparison Plot'));
end

function [input] = plot_margins(input, config)

    for j=1:size(input,2)    
        %index for area of mu_EP_min less margin
        [~, min_index] = min(input(j).result.mu_EP);
        margin_limit =input(j).result.mu_EP_min*(1+config.margin/100);
        
       [~, min_addmargin_index]  = min(abs(input(j).result.mu_EP(min_index:end) - margin_limit ));
       min_addmargin_index= min_addmargin_index+min_index;
            if min_addmargin_index > size(input(j).result.mu_EP ,2)
                min_addmargin_index = size(input(j).result.mu_EP ,2);
            end
        %index for area of mu_EP_min greater margin
       [~, min_lessmargin_index] = min(abs(input(j).result.mu_EP(1:min_index) - margin_limit ));
       
        input(j).mu_min_margin     = input(j).result.mu_EP(min_lessmargin_index:min_addmargin_index);
        input(j).mu_min_margin_c_e = input(j).c_e(min_lessmargin_index:min_addmargin_index);
    end 
end

function [landscape] = make3DFitnesslandscape(input, data, config)
%F_span, c_e_span, m_0, eta_PPU, prop_type,m_struct, I_tot,index, Descriptions)
landscape = struct();

%make x-y grids
for i=1:size(input,2)
    fieldname = strcat('case',num2str(i));
    c_e_min =   config.(char(input(i).thruster_type)).c_e_min;
    c_e_max =   config.(char(input(i).thruster_type)).c_e_max;
    F_min   =   config.(char(input(i).thruster_type)).F_min;
    F_max   =   config.(char(input(i).thruster_type)).F_max;
    
    F_span   = linspace(F_min, F_max, config.res_landscape);  
    c_e_span = linspace(c_e_min, c_e_max,config.res_landscape);
    [landscape .(fieldname).F_grid, landscape .(fieldname).c_e_grid]   = meshgrid(F_span, c_e_span);
end



for i=1:size(input,2)
    
     %   analysis(j).c_e = c_e;
     %   analysis(j).m_dot= input(j).F./c_e;
     %   analysis(j).P_jet_ref = 1/2.*analysis(j).m_dot.*c_e.^2;
     %   analysis(j).PPU_P = analysis(j).P_jet_ref./( analysis(j).PPU_eff*analysis(j).eff);
    
end






figure('Name','Mass Fraction EP System Fitness Landscape');
view(-142,30)
for j=1:index
mu_EP_mesh = [];
    for i=1:size(F_mesh,1)
        F_col = F_mesh(:,i);
        c_e_col = c_e_mesh(:,i);

        mu_EP_col = calculateMassFraction(F_col, c_e_col, m_0, eta_PPU, prop_type,m_struct, I_tot,j);

        mu_EP_mesh = [mu_EP_mesh mu_EP_col];
    end

mesh_handle = mesh(F_mesh,c_e_mesh,mu_EP_mesh, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);%transparent faces, opqaue edges
hold on;
xlabel('F / N');
ylabel('c_e / m/s');
zlabel('\mu_{EP} / - ');

set(mesh_handle, 'FaceColor',map(j,:),'edgecolor',map(j,:));
mesh_collector{end+1} = mu_EP_mesh;

end

legend(Descriptions);

%colormap
map = [ 0.      0.      1.
        0       0.5000  0.
        1.0000  0.      0.
        0       0.7500  0.7500
        0.7500  0       0.7500
        0.7500  0.7500  0
        0.2500  0.2500  0.2500
        0.      1.      0.];
colormap(map);
hold off

%inital min mesh with 1s
mu_min_mesh = ones(size(mesh_collector{1},1),size(mesh_collector{1},2));
index_matrix = ones(size(mesh_collector{1},1),size(mesh_collector{1},2));

%Total minimum mesh
for i=1:size(mesh_collector,2)
    %prop-type
   new_min_mesh = min(mu_min_mesh,mesh_collector{i});
   
    for j=1:size(mesh_collector{1},1)
        for k=1:size(mesh_collector{1},2)
            if new_min_mesh(j,k)<mu_min_mesh(j,k)
                index_matrix(j,k) = i;
            end
        end
    end    
    mu_min_mesh = new_min_mesh;                       %shift reference mesh to get min of all meshes

end


%preallocate
color_matrix_concept = zeros(size(mesh_collector{1},1), size(mesh_collector{1},2),3);
color_matrix_propellant = zeros(size(mesh_collector{1},1), size(mesh_collector{1},2),3);

    for j=1:size(mesh_collector{1},1)
        for k=1:size(mesh_collector{1},2)
            switch index_matrix(j,k)
                case 1 %NH3
                    color_conc = [0. 0. 1.];
                    color_prop = [0. 0.5 0.];
                case 2 %He
                    color_conc = [0    0.5000   0];
                    color_prop = [0.5 0. 0.5];
                case 3 %NH3
                    color_conc = [ 1.0000         0         0];
                    color_prop = [0. 0.5 0.];
                case 4 %He
                    color_conc = [0    0.7500    0.7500];
                    color_prop = [0.5 0. 0.5];
                case 5 %NH3
                    color_conc = [ 0.7500         0    0.7500];
                    color_prop = [0. 0.5 0.];
                case 6 %NH3
                    color_conc = [ 0.7500    0.7500         0];
                    color_prop = [0. 0.5 0.];

                case 7 %HE
                    color_conc = [0.2500    0.2500    0.2500];
                    color_prop = [0.5 0. 0.5];
                case 8 %HE
                    color_conc = [0. 0. 1.];
                    color_prop = [0.5 0. 0.5];
            end
            color_matrix_concept(j,k,:) = color_conc;
            color_matrix_propellant(j,k,:) = color_prop;

        end
    end

figure('Name','Mass Fraction EP System Combined Minimal Fitness Landscape');
mesh(F_mesh,c_e_mesh,mu_min_mesh, color_matrix_concept,'FaceAlpha', 0., 'EdgeAlpha', 1.);
xlabel('F / N');
ylabel('c_e / m/s');
zlabel('\mu_{EP} / - ');
view(-142,30)
legend;

%color according to propellant
figure('Name','Mass Fraction EP System Combined Minimal Fitness Landscape');
mesh(F_mesh,c_e_mesh,mu_min_mesh, color_matrix_propellant,'FaceAlpha', 0., 'EdgeAlpha', 1.);
xlabel('F / N');
ylabel('c_e / m/s');
zlabel('\mu_{EP} / - ');
view(-142,30)

x_mesh = F_mesh;
y_mesh = c_e_mesh;
z_mesh = mu_min_mesh;
colors = color_matrix_propellant;

end

