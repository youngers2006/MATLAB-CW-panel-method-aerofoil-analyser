function [X,Z] = panelgen(NACA_code, N, AoA)
% This function will return two arrays, of X coordinates and Z coordinates respectively. The function accepts inputs of a NACA code formatted as a character string, N the number of pannels to be generated, note that N+1 panels will be generated as an aditional wake panel is also generated, and the angle of attack of the airfoil
    m = str2double(NACA_code(1)) / 100; % the maximum camber (m) is found fro the first character of the NACA code, it must be divided by 100 as it is provided as a percentage 
    p = str2double(NACA_code(2)) / 10; % the location of maximum camber is found by dividing the second character of the NACA code by 10 as it is given in tneths of a chord
    t = str2double(NACA_code(3:4)) / 100; % the maximum aerofoil thickness is found by dividing the last 2 digits of the NACA code by 100 as it is given as a percentage 
    
    % Calculating x array
    i = linspace(1, round((N + 1) / 2), round((N + 1) / 2)); % defines set of points on the x axis
    x = 1 - 0.5 * (1 - cos(2 * pi * (i - 1) / N)); % Using a cosine distribution given in the brief to find position of the panel end points
    
    % Thickness and camber functions
    yt = 5*t*(0.2969*x.^0.5-0.126*x-0.3516*x.^2+0.2843*x.^3-0.1015*x.^4); % calcaulating the thickness function as seen in the brief
    before_max_camber = (x <= p);  % Logical array for x before or at p for use to find camber line
    after_max_camber = ~before_max_camber;  % Logical array for x after p for use to find camber line

    % Preallocate yc and theta to increase computational sppeds by avoiding resizing
    yc = zeros(size(x));
    theta = zeros(size(x));

    % Compute yc and theta using logical indexing
    yc(before_max_camber) = (m / p^2) .* (2 * p * x(before_max_camber) - x(before_max_camber).^2); % calculates yc using the equation in the brief for the camber line, this is the camber line where x <= p
    theta(before_max_camber) = atan((2 * m / p^2) .* (p - x(before_max_camber))); % calculate theta using the equations in the brief where the argument in the tan function is ycprime, these are the theta values where x <= p
    yc(after_max_camber) = (m / (1 - p)^2) .* ((1 - 2 * p) + 2 * p * x(after_max_camber) - x(after_max_camber).^2); % calculates yc using the equation in the brief for the camber line, this is the camber line where x > p
    theta(after_max_camber) = atan((2 * m / (1 - p)^2) .* (p - x(after_max_camber))); % calculate theta using the equations in the brief where the argument in the tan function is ycprime, these are the theta values where x > p

    % Using the equations in the brief the lower and upper coordinate arrays can be found for both x and y using the x values, theta, yc and yt values we have found already
    xU = x - yt .* sin(theta);
    zU = yc + yt .* cos(theta);
    xL = x + yt .* sin(theta);
    zL = yc - yt .* cos(theta);
    
    % force the first point of both x and z arrays through the trailing edge to ensure the airfoil is whole
    xU(1) = 1;
    xL(1) = 1;
    zU(1) = 0;
    zL(1) = 0;
    
    X = [xL,flip(xU,2)]; % concatenates xL and xU, xU is flipped to ensure the data is plotted clockwise
    Z = [zL,flip(zU,2)]; % concatenates zL and zU, zU is flipped to ensure the data is plotted clockwise
    X(2*round((N+1)/2)+1) = 1000000; % creates wake pannel by having an x value at a number that is essentially inifinite when compared to the width of the arifoil
    Z(2*round((N+1)/2)+1) = 1000000 * tan(AoA); % creates wake pannel by having an x value at a number that is essentially inifinite when compared to the width of the arifoil
    
    if N/2 == round(N/2) % when N is even there will be an additional point so this must be removed, the additional point will be at the index N/2+1
        % remove the point at this index by making it empty 
        Z(N/2+1) = []; 
        X(N/2+1) = [];
    end
end 

