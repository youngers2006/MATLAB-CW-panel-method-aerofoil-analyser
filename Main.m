% Housekeeping
clear
clc

% Take in user inputs for the NACA arifoil code, the number of pannels to approximate the airfoil with, the freestream velocity and the angle of attack in degrees
NACA_code_input = input("Input the the NACA code of the aerofoil to be used: ");
N = input("Input number of pannels: ");
Ufs = input("Input the free stream velocity (m/s): ");
AoAdegrees = input("Input the angle of attack (degrees): ");
NACA_code = num2str(NACA_code_input); % converts the NACA_code given as as integer into a string to allow it to be indexed
AoA = (pi/180)*(AoAdegrees); % converts the angle of attack given in degrees to radians

% Calculating Mu, Cl and the discretised airfoil for a standard case
[x,z] = panelgen(NACA_code, N, AoA); % Call panelgen to discretise the airfoil into N+1 pannels (including the wake pannel) and return two arrays, one of the x points and one of the z points to for the discretised airfoil
[Cl,Mu] = Cl_Mu_calculator(x,z,N,AoA,Ufs); % call Cl_Mu_calculator, inputting the coordinates from the panelgen function, the number of pannels, the angle of attack and the freestream velocity. This returns the Cl of the airfoil calculated with these inputs and the N+1 pannel strength values (Mu)
fprintf('Lift Coefficient (C_L): %.4f\n', Cl); % output the Cl value in a formatted string

% calculating flowfield and vector arrows for a standard case
gridpoints_stream = 50; % set the number of points along each axis of the grid to plot the velocity stream, this value must be high to ensure the gird is very fine
gridpoints_vector = 14; % set the number of points along each axis of the grid to plot the velocity vectors, this value can be lower than that used to plot the stream as the grid can be more coarse
[UFlowFieldstream, VFlowFieldstream,Xgridstream, Zgridstream] = velocity_grid(gridpoints_stream,Ufs,AoA,Mu,x,z,N); % Call the velocity grid function to return a grid of horizntal and vertical velocity values and a corresponding grid of coordinates for them to be plotted on for the stream plot
[UFlowFieldvector, VFlowFieldvector,Xgridvector, Zgridvector] = velocity_grid(gridpoints_vector,Ufs,AoA,Mu,x,z,N); % Call the velocity grid function to return a grid of horizntal and vertical velocity values and a corresponding grid of coordinates for them to be plotted on for the velocity vector plot
streamPlotter(x,z,NACA_code,Xgridstream,Zgridstream,UFlowFieldstream,VFlowFieldstream,N,AoAdegrees); % call the stream plotter function to plot a formatted graph containing the velocity stream values calculated above and display the airfoil
vectorPlotter(x,z,NACA_code,Xgridvector,Zgridvector,UFlowFieldvector,VFlowFieldvector,N,AoAdegrees); % call the vector plotter function to plot a formatted graph containing the velocity vector values calculated above and display the airfoil

if NACA_code_input == 2412 % special case for when NACA_code = 2412 as specified on the brief

    % Getting data from XFOIL_DATA text file
    file = fopen("XFOIL_DATA.txt","r"); % open XFOIL_DATA text file in read mode
    textscan(file,'%*[^\n]',7); % skip the first seven lines of text to get to data
    xfoil_values = fscanf(file,'%f',[7,inf])'; % extracts all data from the table in the text file and saves it to xfoil_values
    fclose(file); % close the file when data has been extracted
    
    % Calculating Cl and Mu values
    Narray = [100,200,500]; % define array of N values (number of pannels) to run through the function as specified by the brief
    AoA_array_degrees = [0,1,2,3,4,5,6,7,8,9,10]; % define array of angle of attack values to run through the function as specified by the brief
    AoA_array = (pi/180) * AoA_array_degrees; % convert all angles of attack in the array to radians so they are compataible with the function
    Cl_array = zeros(length(Narray),length(AoA_array)); % Initialise the matrix of Cl values to avoid resizing and imporve computational speed

    for NarrayIndex = 1:length(Narray) % loop through all values of N to be tested
        for AoAIndex = 1:length(AoA_array) % loop through all values of angle of attack to be tested
                [x2412,z2412] = panelgen(NACA_code, Narray(NarrayIndex), AoA_array(AoAIndex)); % as we now have a new angle of attack and new number of pannels the airfoil discitisation must be repeated to ensure that the panel strength calculation is done accurately
                Cl_array(NarrayIndex,AoAIndex) = Cl_Mu_calculator(x2412,z2412,Narray(NarrayIndex),AoA_array(AoAIndex),Ufs); % calculate Cl for the current N, AoA and new calculated panels for the airfoil for every N, AoA in the arrays above
        end 
    end 

    % plots for CL against alpha
    alphaForPlot = linspace(AoA_array(1),AoA_array(end),20); % Create array of finer AoA values to allow for line of best fit plotting
    gradientOfFit = zeros(length(Narray)); % Initialise array to store all gradient values to avoid resizing
    InterceptOfFit = zeros(length(Narray)); % Initialise array to store all intercept values to avoid resizing

    figure % create new figure
    hold on % set hold on
    for i = 1:length(Narray) % loop through all N values in Narray
        if i == 1
            plot(AoA_array_degrees,Cl_array(i,:),'r','LineWidth',1.5) % plot the datapoints from each array of AoA values against Cl with N held constant for this loop
        elseif i == 2
            plot(AoA_array_degrees,Cl_array(i,:),'g','LineWidth',1.5) % plot the datapoints from each array of AoA values against Cl with N held constant for this loop
        elseif i == 3
            plot(AoA_array_degrees,Cl_array(i,:),'b','LineWidth',1.5) % plot the datapoints from each array of AoA values against Cl with N held constant for this loop
        end 

        poly = polyfit(AoA_array_degrees,Cl_array(i,:),1); % create a ploynomial of order one to fit the datapoints just plotted
        Fit = polyval(poly,alphaForPlot); % evaluate the polynomial at finer spaced AoA values
        plot(alphaForPlot,Fit(i)) % plot line of best fit
        gradientOfFit(1,i) = poly(1); % retrieve gradient of the order 1 polynomial (coefficent of the 1st order term)
        InterceptOfFit(1,i) = poly(2); % retrieve intercept of the order 1 polynomial (coefficent of the 0th order term)
    end 
    
    % Plot values from Xfoil
    plot(xfoil_values(:,1),xfoil_values(:,2),'k-','LineWidth',1.5) % plot all values from the data for the Cl against alpha data for Xfoil
    xlabel(sprintf('Angle of attack [%s]',char(176)),"FontSize",15) % Add label to the X axis
    ylabel('Coefficient of lift','FontSize',15) % add label to the z axis
    title('CL vs AoA for a NACA 2412 aerofoil','FontSize',15) % add title to the graph
    xlim([0 10]) % set axis limits to the region we are assessing as specified in the brief
    legend(sprintf('N = 100, Cl = [%i]*alpha + [%i]',gradient(1,1),InterceptOfFit(1,1)),sprintf('N = 200, Cl = [%i]*alpha + [%i]',gradient(1,2),InterceptOfFit(1,2)), sprintf('N = 500, Cl = [%i]*alpha + [%i]',gradient(1,3),InterceptOfFit(1,3)),"Xfoil")
    hold off % set hold off

    %saveas(gcf,sprintf('CL against AoA for NACA %d aerofoil',NACA_code),'fig'); %Save the graph
    
    % Streamline and vector plots for specific case   
    gridpoints_stream2412 = 50; % set the number of points along each axis of the grid to plot the velocity stream, this value must be high to ensure the gird is very fine
    gridpoints_vector2412 = 14; % set the number of points along each axis of the grid to plot the velocity vectors, this value can be lower than that used to plot the stream as the grid can be more coarse
    [x2412,z2412] = panelgen(NACA_code, Narray(end), AoA_array(end)); % Recalculate panel points for end case using panelgen
    [Clend,Mu2412] = Cl_Mu_calculator(x2412,z2412,Narray(end),AoA_array(end),Ufs); % Recalculate panel strength values using Cl_Mu_calculator
    [UFlowFieldstream2412, VFlowFieldstream2412,Xgridstream2412, Zgridstream2412] = velocity_grid(gridpoints_stream2412,Ufs,AoA_array(end),Mu2412,x,z,Narray(end)); % Use velocity grid function to calculate velocity grids for streamlines
    [UFlowFieldvector2412, VFlowFieldvector2412,Xgridvector2412, Zgridvector2412] = velocity_grid(gridpoints_vector2412,Ufs,AoA_array(end),Mu2412,x,z,Narray(end)); % Use velocity grid function to calculate velocity grids for velocity vectors
    streamPlotter(x,z,NACA_code_input,Xgridstream2412,Zgridstream2412,UFlowFieldstream2412,VFlowFieldstream2412,N,AoAdegrees); % Use stream plotter to plot streamlines for the end case
    vectorPlotter(x,z,NACA_code_input,Xgridvector2412,Zgridvector2412,UFlowFieldvector2412,VFlowFieldvector2412,N,AoAdegrees); % use vector plotting function to plot vectors for end case
end 

% Functions used in the code - I chose to make these functions as they are functions that are repeated several times throughout the script so to improve code readability and reduce the amount of code in the script and reduce the possibility of errors I have used these Functions
function [Cl,Mu] = Cl_Mu_calculator(x,z,N,AoA,Ufs)
% Function to calculate the Coefficient of lift and the panel strengths for a given set of panel coordinates, panel number, angle of attack and freestream velocity
    % Initialise A and B vectors to solve for the panel strengths
    A = zeros(N+1, N+1); % A will be an N+1 by N+1 matrix as there is N+1 panels including the wake panel that will all influence eachother so there should be N+1 columns and there should be N+1 rows as N panels are solved for and the wake panel is solved for via the Kutta condition
    B = zeros(N+1, 1); % B will be a vector of length N+1 as there is N+1 panels considered so there should be N+1 solutions for the N+1 equations in the system
    for i = 1:N % loop through all pannel midpoints (our regions of interest)
        beta = atan2(z(i+1) - z(i), x(i+1) - x(i)); % calculate the pannel angle beta using the equation for beta provided in the brief
        for j = 1:N+1 % loop through all pannel bounds, as we are considering superposition all pannel strengths effect all pannels
            p = [(x(i+1)+x(i))/2,(z(i+1)+z(i))/2]; % coordinates of the ith pannels midpoint
            p1 = [x(j),z(j)]; % coordinates of the start point of the ith pannel
            p2 = [x(j+1),z(j+1)]; % coordinates of the endpoint of the ith pannel
            [u, v] = cdoublet(p, p1, p2); % calculate the induced velocities u and v as a results of the strengths of the current pannels of interest 
            A(i, j) = (v * cos(beta)) - (u * sin(beta)); % using the equation from the coursework brief this is one entry in the A matrix which represents a system of equations
        end 
        B(i,1) = -1 * Ufs * sin(AoA - beta); % calculate the ith entry for B the RHS of the equation 
    end 

    % applying the kutta condition for the wake panel to allow the system of equations to be solved
    A(N+1,N+1) = 1;
    A(N+1,1) = 1;
    A(N+1,N) = -1;
    B(N+1) = 0;

    Mu = A\B; % solve the system of equations to now have a vector of pannel strengths Mu for all pannels 
    Cl = -2 * Mu(end) / Ufs; % calculate the lift coefficient of the airfoil using th equation given in the brief
end 

function [UFlowField, VFlowField,Xgrid, Zgrid] = velocity_grid(gridpoints,Ufs,AoA,Mu,x,z,N)
    % calculating flowfield and vector arrows
    [Xgrid, Zgrid] = meshgrid(linspace(-0.2, 1.2, gridpoints), linspace(-0.7, 0.7, gridpoints)); % use meshgrid to create a matrix of points (a grid) for x and z with the coarseness of the grid dependent on the number of points specified on each axis
    % Initialize velocity fields to improve computing time by avoiding resizing
    UFlowField = zeros(size(Xgrid));
    VFlowField = zeros(size(Zgrid));

    % Compute velocity field
    for i = 1:numel(Xgrid) % loops from one to the number of elements in the matrix Xgrid
        Uflow = Ufs * cos(AoA); % calculate first term in the velocity equation in the brief for horizontal velocity
        Vflow = Ufs * sin(AoA); % calculate first term in the velocity equation in the brief for vertical velocity
        p = [Xgrid(i), Zgrid(i)]; % Grid point of interest
        for k = 1:N+1 % loop through all panels
            [ugrid, vgrid] = cdoublet(p, [x(k), z(k)], [x(k+1), z(k+1)]); % calculate the induced velocity at the grid point of interest by each individual panel with cdoublet function
            Uflow = Uflow + Mu(k) * ugrid; % calculate the effeect of the kth panel on the horizontal velocity at the grid point of interest
            Vflow = Vflow + Mu(k) * vgrid; % calculate the effeect of the kth panel on the vertical velocity at the grid point of interest
        end
        UFlowField(i) = Uflow; % add the value of horizonal velocity at the point of interest to the horizontal velocity value grid
        VFlowField(i) = Vflow; % add the value of vertical velocity at the point of interest to the vertical velocity value grid
    end

    % Mask points inside the airfoil
    inpoly = inpolygon(Xgrid, Zgrid, x, z); % checks all gridpoints against airfoil plot points to check if they lie within the airfoil
    UFlowField(inpoly) = NaN; % if it is the case that they lie within the airfoil they are set to value NaN and so are not plotted
    VFlowField(inpoly) = NaN; % if it is the case that they lie within the airfoil they are set to value NaN and so are not plotted
end 

function streamPlotter(x,z,NACA_code,Xgrid,Zgrid,UFlowField, VFlowField,N,AoAdegrees)
% This function is to plot the velocity streamlines of the flow for a given airfoil
    figure % create new figure to avoid overlap
    hold on % set hold on
    plot(x(1:end-1),z(1:end-1),'r-','Linewidth',2) % plot arifoil, excluding the wake panel
    ylim([-0.7 0.7]) % limit x axis to keep the streamlines inside
    xlim([-0.2 1.2]) % limit y axis to keep the streamlines inside
    title(sprintf('Streamline field for NACA %i at AoA = %f degrees, and N = %i',NACA_code, AoAdegrees, N)); % Give title to plot
    xlabel('x'); % label x axis
    ylabel('z'); % label y axis
    streamlines = streamslice(Xgrid, Zgrid, UFlowField, VFlowField); % plot streamlines
    set(streamlines,'color','b','linewidth',1); % format streamlines
    grid on; % set a grid
    hold off
    saveas(gcf,sprintf('Streamlines plot for NACA %d aerofoil',NACA_code),'fig'); % Save the graph
end 

function vectorPlotter(x,z,NACA_code,Xgrid, Zgrid, UFlowField, VFlowField,N,AoAdegrees)
% This function is to plot the velocity vecotrs of the flow for a given airfoil
    figure % create new figure to avoid overlap
    hold on % set hold on
    plot(x(1:end-1),z(1:end-1),'r-','Linewidth',2) % plot arifoil, excluding the wake panel
    ylim([-0.7 0.7]) % limit x axis to keep the streamlines inside
    xlim([-0.2 1.2]) % limit y axis to keep the streamlines inside
    title(sprintf('vector field for NACA %d at AoA = %d degrees, and N = %d',NACA_code, AoAdegrees, N)); % Give title to plot
    xlabel('x'); % label x axis
    ylabel('z'); % label y axis
    quiver(Xgrid, Zgrid, UFlowField, VFlowField,'Color','b',LineWidth=1); % plot all velocity vectors
    grid on; % set a grid
    hold off
    saveas(gcf,sprintf('velocity vector plot for NACA %d aerofoil',NACA_code),'fig'); %Save the graph
end 