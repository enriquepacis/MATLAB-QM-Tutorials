close all;
clear all;

%% Constants
hbar = 6.58E-16; % [eV*s] the reduced Planck constant

%% Parameters
N = 75; % [] number of quantum dots
nt = 201; % [] number of time steps

gamma = 0.05; % [eV] transition (tunneling) energy


%% Calculations

% Make the Hamiltonian in the site basis

H = diag( -gamma * ones(1, N-1), +1 ) + diag( -gamma * ones(1, N-1), -1 );


% Calculate the eigenstates of H (in the site basis)

% E_n: ordered eigenvalues of the system
% phi_n: ordered eigenstates of the system
[phi_n, E_n] = eig(H); 

Selection = 4;

phi_sel = phi_n(:, Selection);

prob_dist = conj(phi_sel) .* phi_sel;


% Time evolution
% psi_0 = phi_n(:,1); % select the ground state as the initial state
psi_0 = (1/sqrt(2)) * (phi_n(:,2) + phi_n(:,3)); % select the ground state as the initial state

Eo = E_n(2,2) - E_n(1,1); % [eV] energy scale for system
To = hbar/Eo; % [s] a time scale for the system

t = linspace(0, 2*To, nt);

psi_t = zeros(N, nt);

for t_idx = 1:nt % Calculate psi_t at each time step
    Ut = expm( - (1i * t(t_idx) / hbar) * H );
    psi_t(:, t_idx) = Ut * psi_0;
end
    
prob_dist_t = conj(psi_t).*psi_t;


%% Visualization

bar(prob_dist);
set(gca, 'FontName', 'Times', 'FontSize', 20);
xlabel('Site')
ylabel('Probability')
grid on;

figure;
surf(ones(N,1) * (1/To) * t , (1:N)' * ones(1,nt) , prob_dist_t);
shading interp
set(gca, 'FontName', 'Times', 'FontSize', 20);
ylabel('Site');
xlabel('$t/T_o$', 'Interpreter', 'latex');
zlabel('Probability');


