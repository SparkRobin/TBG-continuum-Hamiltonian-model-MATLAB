% Twisted Bilayer Graphene Band Structure (Fixed Truncation)
% Fig 3a of PNAS 2011:
% <Moiré bands in twisted double-layer graphene> by Rafi Bistritzer and Allan H. MacDonald

%% 1. Parameters and Constants (User's Original Settings)
close all
clear
tic

% if consider lattice relaxation, it is appropriate to adopt wAA/wAB=0.8
% from <PHYS. REV. X 8, 031087 (2018)>

w_AB = 0.110; % eV
w_AA = 0.110;
% w_AA=w_AB*0.8;

% Magic Angle: 1.05deg
theta_deg = 1.05;
theta = deg2rad(theta_deg);

a0 = 2.46;                  % Lattice constant (Angstroms)
d = 1.42;
t = 2.97;
hvF = 3/2*t*d;              % Fermi Velocity

% Truncation Parameter
% The smaller the twist angle, the bigger the tructions
if theta_deg > 2
    truncation_shell = 3;
elseif theta_deg > 0.8
    truncation_shell = 10;  % For magic angle ~1.05
else
    truncation_shell = 20; % For 0.5 deg (需消耗更多内存和计算时间)
end

%% 2. Geometry and Vectors

% Monolayer K point magnitude
K_mag = 4 * pi / (3 * a0);

% Moiré modulation wavevector
kD = 2 * K_mag * sin(theta / 2);

% Coupling Vectors q (shifts from L1 to L2)
qb  = kD * [0; -1];
qtr = kD * [sqrt(3)/2; 1/2];
qtl = kD * [-sqrt(3)/2; 1/2];

% Moiré Reciprocal Basis Vectors
bM1 = sqrt(3) * kD * [1/2; -sqrt(3)/2];
bM2 = sqrt(3) * kD * [1/2; sqrt(3)/2];

%% 3. Dynamic Grid Generation

basis_sites = []; % [n, m, Layer_Index]
count = 0;
scan_range = -truncation_shell:truncation_shell;

for n = scan_range
    for m = scan_range
        % the momentum-space lattice is truncated at hexagonal region
        if abs(n) + abs(m) + abs(n-m) <= 2*truncation_shell

            % Layer 1 Site
            count = count + 1;
            basis_sites(count, :) = [n, m, 1];

            % Layer 2 Site
            count = count + 1;
            basis_sites(count, :) = [n, m, 2];
        end
    end
end

num_sites = size(basis_sites, 1);
disp(['Total Matrix Size: ', num2str(num_sites * 2)]);

% Compute Q shift vectors relative for K point of layer 1
Q_vecs = zeros(2, num_sites);
for i = 1:num_sites
    n = basis_sites(i, 1);
    m = basis_sites(i, 2);
    layer = basis_sites(i, 3);

    %  G = n*bM1 + m*bM2
    vec = n * bM1 + m * bM2;

    if layer == 2 %Momentum Shift
        vec = vec + qb;
    end
    Q_vecs(:, i) = vec;
end

%% 4. Coupling Matrices
phi = 2*pi/3;

% T matrices 
T1 = [w_AA, w_AB;
    w_AB, w_AA];
T2 = [w_AA * exp(1i*phi), w_AB;
    w_AB * exp(-1i*phi), w_AA * exp(1i*phi)];
T3 = [w_AA * exp(-1i*phi), w_AB;
    w_AB * exp(1i*phi), w_AA * exp(-1i*phi)];

%% 5. High symmetry k-path Construction 
Cpoint = bM2;
Apoint = (2*bM1 + bM2) / 3;
Mm = (bM1 + 2*bM2) / 2;
Bpoint = (bM1 + 2*bM2) / 3;
Dpoint = (bM1 + bM2);

% Path: A->B->C->D->M->A
num_pts = 20;
path1 = generate_path(Apoint, Bpoint, num_pts);
path2 = generate_path(Bpoint, Cpoint, num_pts);
path3 = generate_path(Cpoint, Mm, num_pts);
path4 = generate_path(Mm, Dpoint, num_pts);
path5 = generate_path(Dpoint, Apoint, num_pts);

k_vecs = [path1, path2, path3, path4, path5];
k_dist = [0, cumsum(sqrt(sum(diff(k_vecs, 1, 2).^2)))];

%% 6. Pre-calculate Couplings

couplings = []; % [index_i, index_j, type]
tolerance = 1e-5;

fprintf('Building coupling map...\n');
for i = 1:num_sites
    if basis_sites(i, 3) == 1 % Source is Layer 1
        for j = 1:num_sites
            if basis_sites(j, 3) == 2 % Target is Layer 2
                dq = Q_vecs(:, j) - Q_vecs(:, i); % momentum difference

                % if momentum difference equals qj
                if norm(dq - qb) < tolerance
                    couplings = [couplings; i, j, 1]; % Type 1 (T1)
                elseif norm(dq - qtr) < tolerance
                    couplings = [couplings; i, j, 2]; % Type 2 (T2)
                elseif norm(dq - qtl) < tolerance
                    couplings = [couplings; i, j, 3]; % Type 3 (T3)
                end
            end
        end
    end
end

%% 7. Hamiltonian Loop
bands = zeros(num_sites*2, size(k_vecs, 2));

% Rotated Pauli Matrices
sx = [0 1; 1 0];
sy = [0 -1i; 1i 0];
theta1 = -theta/2;
theta2 = theta/2;
sx1 = sx*cos(theta1) - sy*sin(theta1);
sy1 = sx*sin(theta1) + sy*cos(theta1);
sx2 = sx*cos(theta2) - sy*sin(theta2);
sy2 = sx*sin(theta2) + sy*cos(theta2);

fprintf('Calculating bands...\n');
for idx = 1:size(k_vecs, 2)
    k_curr = k_vecs(:, idx);

    %  q_base = 0 when k_curr = Apoint.
    q_base = k_curr - Apoint;

    H = zeros(num_sites*2, num_sites*2);

    % 1. Diagonal (Kinetic)
    for n = 1:num_sites
        q_loc = q_base + Q_vecs(:, n);

        if basis_sites(n, 3) == 1
            h_block = hvF * (q_loc(1) * sx1 + q_loc(2) * sy1);
        else
            h_block = hvF * (q_loc(1) * sx2 + q_loc(2) * sy2);
        end

        r_idx = (n-1)*2+1;
        H(r_idx:r_idx+1, r_idx:r_idx+1) = h_block;
    end

    % 2. Off-Diagonal (Tunneling) - Using pre-calculated map
    for c = 1:size(couplings, 1)
        i = couplings(c, 1);
        j = couplings(c, 2);
        type = couplings(c, 3);

        switch type
            case 1, T = T1;
            case 2, T = T2;
            case 3, T = T3;
        end

        r = (i-1)*2+1;
        c_idx = (j-1)*2+1;

        H(r:r+1, c_idx:c_idx+1) = T;
        H(c_idx:c_idx+1, r:r+1) = T'; % Hermitian
    end

    bands(:, idx) = sort(real(eig(H)));
    fprintf('Calculating Eigenvalues for k-point: %d / %d \n', idx, size(k_vecs, 2));

end

fprintf('Band structure calculated...\n');
%% 8. Plotting
figure('Color', 'w');

% Select bands near energy 0 (Charge Neutrality)
mid_band = num_sites;
num_plot = 14; % Plot 14 bands closest to the Dirac point
range_idx = (mid_band - num_plot/2 + 1) : (mid_band + num_plot/2);
band_plot=bands(range_idx, :);
plot(k_dist, band_plot * 1000, '-', 'LineWidth', 1.2);

xlim([0, max(k_dist)]);

ylabel('E (meV)');
xlabel('Wave Vector');

% X-ticks
tick_indices = [1, num_pts+1, 2*num_pts+1, 4*num_pts+1, 5*num_pts];
xticks(k_dist(tick_indices));
xline(k_dist(tick_indices), '--', 'Color', [0.7 0.7 0.7]);
xticklabels({'A', 'B', 'C', 'M', 'D', 'A'});

title(['TBG Band Structure (\theta = ' num2str(theta_deg) '^\circ)']);
box on;

filename_png = sprintf('Eband_theta%.2fdeg.png',theta_deg);
saveas(gcf,filename_png , 'png');

filename_mat = sprintf('Eband_theta%.2fdeg.mat',theta_deg);
save(filename_mat);

toc
%% Functions
function pts = generate_path(p1, p2, n)
pts = p1 + (p2 - p1) * linspace(0, 1, n);
end