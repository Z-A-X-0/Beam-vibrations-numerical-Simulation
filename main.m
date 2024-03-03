%% Instructions: %%
%%% General: Remove the % in front of the functions in lines 64-103 to
%   output a function %%%
% 1. Please output Project Part 1 and Part 2 separately
% 2. In the first part of the project, the functions must be output in sequence
%    for the results of the functions to be adopted, otherwise, errors will occur
% 3. Please do not output tasks that involve plotting at the same time, as
%    errors will occur in the plot
% 4. Please output the functions of Task 14 (Eges) individually, otherwise,
%    errors will occur in the plots

%% Initializing all values within the project
n=3; % Number of elements in 1
hat_n=1; % Number of additional evaluation points per element in 1
tilde_n=7; % Order of quadrature in 1
nP=100; % Number of time steps in 1
beta=1/4; % Newmark coefficient in 1
gamma=1/2; % Newmark coefficient in 1
eta=0.1; % Time step size in s
l=1; % Length of the beam in m
mu=1; % Linear mass density in kg/m
E=1; % Elastic modulus in N/m^2
I=1; % Moment of inertia in m^4
q=1; % Distributed load in N/m
B=[

0 1 0; % Displacement at the left end in m
0 2 0; % Slope at the left end in 1
n 3 0; % Moment at the right end in Nm
n 4 0 % Shear force at the right end in N
];

h = l/n;  

params = {
    0.25, 0.5, 0.1; % {beta, gamma, eta}
    0.25, 0.5, 1;
    0.5, 1, 0.1;
    0.5, 1, 1
};

% % Task 15
% 
% % Linear mass density in kg/m
% mu = @(x) 1;
% 
% % Elastic modulus in N/m^2
% E = @(x) 1;
% 
% % Moment of inertia in m^4
% I = @(x) 1;
% 
% % Distributed load in N/m
% q = @(x) 1;


%% List of all functions:
% To output % remove % for function
%% First project part
%getindizes(n)
%Mbar = getMbar(h, mu,n)
%Sbar = getSbar(E, I, h, n)
%qbar = getqbar(q, h, n)
%M = getM(h, mu,n)
%S = getS(E, I, h, n)
%vq = getvq(q, h, n)
%C = getC(B,n)
%vn = getvn(B,n)
%vd = getvd(B)
%Me = getMe(h, mu, n)
%Se = getSe(E, I, h, B, n)
%ve = getve(q, h, B, n)
%aE = getaE(E, I, h, B, q, n)
%getplot(aE, n, l)
%solveNewmark(nP, E, I, h, B, n, l, mu, q, eta, beta, gamma)
%getrelFehler(mu, q, E, I, l)

%%% Please output each individually, otherwise errors will occur in the plot %%%
%%% Additionally, please either: comment out the plot() and pause() command in the solveNewmark function (line 402 and 403),
%%% or: close the plots of the bending lines (might take some time), the last plot is then the desired one %%%
%Eges(n, B, q, h, E, I, mu, nP, l, params{1,1}, params{1,2}, params{1,3});
%Eges(n, B, q, h, E, I, mu, nP, l, params{2,1}, params{2,2}, params{2,3});
%Eges(n, B, q, h, E, I, mu, nP, l, params{3,1}, params{3,2}, params{3,3});
%Eges(n, B, q, h, E, I, mu, nP, l, params{4,1}, params{4,2}, params{4,3});
%% Second project part
%s = getstencil(tilde_n)
%getphi(tilde_n)
%getddphi(tilde_n)
%gethl(n, l)
%getTinV(n, l, tilde_n)
%getexp(n)
% Aufgaben 19 und 20 funktionieren nicht sind aber geschrieben worden
%% End of function list
%% First part of the project

function [i ,j ,a2l_j, a2l_i, ib, lb, b2l_i] = getindizes(n)


% Laufindizes

L = 0:n-1; % l ∈ L = {0, ... , n-1}
I = 0:3; % i, j ∈ I = {0, 1, 2, 3}

% 3D-Arrays
[i, j, l] = ndgrid(I,I,L);


a2l_i = 2 * l + i;
a2l_j = 2 * l + j;

% Ausgabe 3D
% disp("3D-Arrays:");
% 
% disp("l:");
% disp(l);
% 
% disp("i:");
% disp(i);
% 
% disp("j:");
% disp(j);
% 
% disp("2l + i:");
% disp(a2l_i);
% 
% disp("2l + j:");
% disp(a2l_j);


% 2D-Arrays
[ib,egal,lb] = ndgrid(I, 1, L);

b2l_i = 2 * lb + ib;

% Ausgabe 2D
%disp("2D-Arrays:");

%disp("l:");
%disp(lb);

%disp("i:");
%disp(ib);

%disp("2l + i:");
%disp(b2l_i);


end

function M_bar = getMbar(h, mu, n)

% Element matrix M
M1 = (mu * h / 420) * [156 22*h 54 -13*h;
22*h 4*h^2 13*h -3*h^2;
54 13*h 156 -22*h;
-13*h -3*h^2 -22*h 4*h^2];

% Erstellen des 3D-Arrays
M_bar = repmat(M1, [1, 1, n]);
end

function S_bar = getSbar(E, I, h, n)

% Element matrix S
S1 = (E * I / h^3) * [12 6*h -12 6*h;
6*h 4*h^2 -6*h 2*h^2;
-12 -6*h 12 -6*h;
6*h 2*h^2 -6*h 4*h^2];

% 3D-Arrays
S_bar = repmat(S1, [1, 1, n]);
end

function q_bar = getqbar(q, h, n)

%Element vektor q
q1 = (q*h / 12) * [6; h; 6; -h];

% 3D-Arrays

q_bar = repmat(q1, [1, 1, n]);
end

function M = getM(h, mu, n)
    M1 = getMbar(h, mu, n); % Assuming getMbar returns a 3D array
    [i, j, a2l_j, a2l_i, ~, ~, ~] = getindizes(n);

    i = a2l_i + 1;
    j = a2l_j + 1;

    m = max(i(:));
    N = max(j(:));

    % Linearizing the indices and values
    lin_i = i(:);
    lin_j = j(:);
    lin_M1 = M1(:);

    % Creating a single sparse matrix
    M = sparse(lin_i, lin_j, lin_M1, m, N);
%     M = full(M); % Converting to a full matrix
end

function S = getS(E, I, h, n)
    S1 = getSbar(E, I, h, n); % Assuming getSbar returns a 3D array
    [i, j, a2l_j, a2l_i, ~, ~, ~] = getindizes(n);

    i = a2l_i + 1;
    j = a2l_j + 1;

    M = max(i(:));
    N = max(j(:));

    % Linearizing the indices and values
    lin_i = i(:);
    lin_j = j(:);
    lin_S1 = S1(:);

    % Creating a single sparse matrix
    S = sparse(lin_i, lin_j, lin_S1, M, N);
%     S = full(S); % Converting to a full matrix
end

function vq = getvq(q, h, n)
    q1 = getqbar(q, h, n); % Assuming getqbar returns a 3D array
    
    [~ ,~, ~, ~, ib, ~, b2l_i] = getindizes(n);
    j = 1;
    i = b2l_i + 1;
    M = max(i(:));
    N = 1;
    
    % Linearizing the indices and values
    lin_i = i(:);
    lin_j = repmat(j, numel(lin_i), 1); % j is 1, so we replicate it to match the length of lin_i
    lin_q1 = q1(:);
    
    % Creating a single sparse matrix
    vq = sparse(lin_i, lin_j, lin_q1, M, N);
end

function C = getC(B,n)

[i ,j ,a2l_j, a2l_i, ib, lb, b2l_i] = getindizes(n);
M = 2 * n +2;  %2*(n+1)

V1 = sum(B(:,2) == 1);
N1 = V1;

xk1 = B(1:V1,1);
E1 = sparse(2*xk1 + 1,i(1:V1)+1,1,M,N1);


V2 = sum(B(:,2) == 2);
N2 = V2;

xk2 = B(V1+1:V1+V2,1);
E2 = sparse(2*xk2 + 2,i(1:V2)+1,1,M,N2);

C = [E1 , E2];

end

function vn = getvn(B, n)
    % Extracting rows from B where the second column is 3 or 4
    rows3 = B(B(:, 2) == 3, :);
    rows4 = B(B(:, 2) == 4, :);
    
    % Creating sparse matrices E3 and E4
    E3 = sparse(2 * rows3(:, 1) + 2, 1:size(rows3, 1), 1, 2 * n + 2, size(rows3, 1));
    E4 = sparse(2 * rows4(:, 1) + 1, 1:size(rows4, 1), 1, 2 * n + 2, size(rows4, 1));
        
   
    % Extracting c3 and c4 vectors from the third column of rows3 and rows4
    c3 = rows3(:, 3);
    c4 = rows4(:, 3);
    
    % Calculating vn
    vn = E3 * c3 + E4 * c4;
end

function vd = getvd(B)

V1 = sum(B(:,2) == 1);

V2 = sum(B(:,2) == 2);

a = B(1:V1,3);

b = B(V1+1:V1+V2,3);

vd = [a ; b];

end

function Me = getMe(h, mu, n)

M = getM(h, mu,n);
Me12 = sparse(2*n +2,2);
Me21 = sparse(2,2*n +2);
Me22 = sparse(2,2);

Me = [M, Me12; Me21, Me22];
end

function Se = getSe(E, I, h, B, n)

S = getS(E, I, h, n);
C = getC(B,n);
Ct = transpose(C);
Se22 = sparse(2,2);

Se = [S, C; Ct, Se22];
end

function ve = getve(q, h, B, n)
vq = getvq(q, h, n);
vn = getvn(B,n);
vd = getvd(B);

ve11 = vq + vn;
ve = [ve11; vd];
end

function aE = getaE(E, I, h, B, q, n)

Se = getSe(E, I, h, B, n);
ve = getve(q, h, B, n);

aE = Se \ ve;

end

function getplot(aE, n, l)

    N = 2 * n + 2;
    a = aE(1:N);
    a2k = a(1:2:end);
    x = linspace(0, 1, n + 1);

    plot(x, a2k)
    title('Biegelinie');
    xlabel('x in m');
    ylabel('w in m');
    yticks();
    ylim([-0.25, 0.25]);
    xlim([0, l]);
    drawnow % Force MATLAB to update the figure immediately
end

function [vE, vE_dot] = solveNewmark(nP, E, I, h, B, n, l, mu, q, eta, beta, gamma)
    % Precompute matrices and vectors
    Se = getSe(E, I, h, B, n);
    Me = getMe(h, mu, n);
    ve = getve(q, h, B, n);

    % Initial conditions
    aE = Se \ ve;
    vE0 = aE;
    vE0_dot = zeros(size(vE0));
    vE0_ddot = zeros(size(vE0));

    % Initialization
    q = 0;
    ve = getve(q, h, B, n);

    % Arrays to store results
    vE = zeros(length(vE0), nP);
    vE_dot = zeros(length(vE0_dot), nP);
    vE_ddot = zeros(length(vE0_ddot), nP);

    % Store initial conditions
    vE(:, 1) = vE0;
    vE_dot(:, 1) = vE0_dot;
    vE_ddot(:, 1) = vE0_ddot;

    % Time-stepping loop
    for p = 1:nP
        vEp_star = vE(:, p) + eta * vE_dot(:, p) + (1/2 - beta) * eta^2 * vE_ddot(:, p);
        vEp_dot_star = vE_dot(:, p) + (1 - gamma) * eta * vE_ddot(:, p);
        
        A = Me + Se * beta * eta^2;
        b = ve - Se * vEp_star;
        vEp1_ddot = A \ b;

        vE(:, p+1) = vEp_star + beta * eta^2 * vEp1_ddot;
        vE_dot(:, p+1) = vEp_dot_star + gamma * eta * vEp1_ddot;
        vE_ddot(:, p+1) = vEp1_ddot;

        % Plot results
        getplot(vE(:, p+1), n, l);
        pause(0.1);  % Pause to create animation effect
    end
end

function getrelFehler(mu, q, E, I, l)
    n_values = [1,2,3,5,7,9,10,15,20,25,30,35,40,45,50,60,70,80,90,100,125,150,175,200,225,250,275,300,500,700,1000];
    error_L2_values = zeros(size(n_values));
    
    % Instead of for-loop, use arrayfun to compute error for each n value
    error_L2_values = arrayfun(@(n) computeError(n, mu, q, E, I, l), n_values);
        
  % Plotting the relative error in a linear diagram
    figure;
    plot(n_values, error_L2_values, 'b-');
    xlabel('n');
    ylabel('error\_L2');
    title('Relative Error over the Node Number');
    grid on;
        
    % Plotting the relative error in a log-log diagram
    figure;
    plot(log(n_values), log(error_L2_values), 'r-');
    xlabel('log n');
    ylabel('log error\_L2');
    title('Relative Error over the Node Number in log-log Diagram');
    grid on;
end
    
function error = computeError(n, mu, q, E, I, l)
    B=[0 1 0; 0 2 0; n 3 0; n 4 0];
    h = l/n;
    M = getM(h, mu, n);
    A = M;
    %   [i ,j ,a2l_j, a2l_i, ib, lb, b2l_i] = getindizes(n); % Here is where getindizes is used.
    
    x = transpose(linspace(0, l, 2*n+2));
    w_x = q / (E * I) * (x.^4 / 24 - l * x.^3 / 6 + l^2 * x.^2 / 4);
    
    aE = getaE(E, I, h, B, q, n);
    alpha = aE(1:2*n+2); % Extracting alpha values
    
    error = sqrt((w_x - alpha)' * A * (w_x - alpha)) / (sqrt(w_x' * A * w_x));
end

function Eges(n, B, q, h, E, I, mu, nP, l, beta, gamma, eta)
    [vE, vE_dot] = solveNewmark(nP, E, I, h, B, n, l, mu, q, eta, beta, gamma);
    
    q = 0;
    % Initialize other matrices and vectors
    C = getC(B, n);
    vn = getvn(B, n);
    vq = getvq(q, h, n);
    Se = getSe(E, I, h, B, n);
    Me = getMe(h, mu, n);
    
    % Calculate Eges without using loops
    calculateSingleEges = @(vE, vE_dot) 0.5 * vE_dot' * Me * vE_dot ...
                                       + 0.5 * (Se(1:2*n+2, :) * vE - vq - C * sparse(2,1) - vn)' * vE(1:2*n+2, :);
    Eges = arrayfun(@(p) calculateSingleEges(vE(:, p), vE_dot(:, p)), 1:nP+1);
    
    % Plot Eges
    timeVector = linspace(0, nP, length(Eges)); % Corrected timeVector calculation
    figure;
    plot(timeVector, Eges);
    xlabel('Time (s)');
    ylabel('Gesamt Energie in (J)');
    xlim([-5 100]); % Set the limits for the y-axis
    ylim([-1e-3 max(Eges)+1e-3]); % Set the limits for the y-axis
    title(['Gestamt Energie ueber die Zeit, \beta=', num2str(beta), ', \gamma=', num2str(gamma), ', \eta=', num2str(eta)]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Second part of the project
function s = getstencil(tilde_n)
    % Generate the support points (x-coordinates) in the interval [0, 1].
    x_stuetzstellen = linspace(0, 1, tilde_n);
    
    % Use ndgrid to create the coordinates for the support points.
    [x, y] = ndgrid(x_stuetzstellen, 0);
    
    % Create the Vandermonde matrix from the support point coordinates.
    V = x.^((0:tilde_n - 1))
    
    i = (0:tilde_n - 1)';
    Vplus1 = 1./(i + 1)
    s = (transpose(V) \ Vplus1)';
    
    s_3D = reshape(s, [1, 1, tilde_n]) % Reshape to 3D array
    s_4D = reshape(s, [1, 1, 1, tilde_n])  % Reshape to 4D array
end

function getphi(tilde_n)
    % Values of the support points
    xk = linspace(0, 1, tilde_n+1);
    
    phi0 = 1 - 3*xk.^2 + 2*xk.^3;
    phi1 = xk - 2*xk.^2 + xk.^3;
    phi2 = 3*xk.^2 - 2*xk.^3;
    phi3 = -xk.^2 + xk.^3;
    
    phi =[phi0; phi1; phi2; phi3]
    phi_4D = reshape(phi, [1, 1, 4, tilde_n+1])  % Reshape to 4D array
end

function getddphi(tilde_n)
    % Values of the support points
    xk = linspace(0, 1, tilde_n+1);
    
    ddphi0 = -6 + 12*xk;
    ddphi1 = -4 + 6*xk;
    ddphi2 = 6 -12*xk;
    ddphi3 = -2 + 6*xk;
    
    ddphi =[ddphi0; ddphi1; ddphi2; ddphi3]
    ddphi_4D = reshape(ddphi, [1, 1, 4, tilde_n+1])  % Reshape to 4D array
end

function hl = gethl(n, l)
    xl_plus_1 = l + l/n;
    xl = l;
    
    hl = ones(n,1) .* (xl_plus_1-xl);
    hl_3D = reshape(hl, [1, 1, n])  % Reshape to 3D array
end

function getTinV(n, l, tilde_n)
    hl = gethl(n, l);
    
    xl = linspace(0, 2*l/n, n);
    xk = linspace(0, 1, tilde_n+1);
    
    Tinv = hl * xk + xl'
    Tinv_4D = reshape(Tinv, [n, 1, 1, tilde_n+1])  % Reshape to 4D array
end

function getexp(n)
[i ,j ,a2l_j, a2l_i, ib, lb, b2l_i] = getindizes(n);

 
    DreiD_Array = [(i(1:2,1:2,:) + j(1:2,1:2,:)),(i(1:2,1:2,:) + j(1:2,1:2,:));...
    (i(1:2,1:2,:) + j(1:2,1:2,:)),(i(1:2,1:2,:) + j(1:2,1:2,:))]

    ZweiD_Array = DreiD_Array(:,1,:)
end

