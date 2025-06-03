close all; clear; clc

figure; hold on; grid on

%% draw GMRES results
evals = readmatrix('GMRES_evals.txt');
residuals = readmatrix('GMRES_residual.txt');
plot(evals(:,end), residuals(:,end), 'LineWidth', 2)

%% draw IDR(s) results
evals = readmatrix('IDRs_evals.txt');
residuals = readmatrix('IDRs_residual.txt');
plot(evals(:,end), residuals(:,end), 'LineWidth', 2)

%% draw BiCGSTAB results
evals = readmatrix('BiCGSTAB_evals.txt');
residuals = readmatrix('BiCGSTAB_residual.txt');
plot(evals(:,end), residuals(:,end), 'LineWidth', 2)

%% plot annotations
yline(1e-6)
set(gca, 'YScale', 'log')
xlabel('matrix-vector operations', 'FontSize', 20, 'Interpreter', 'latex')
ylabel('$\|A\mathbf{x}_n-\mathbf{b}\|/\|\mathbf{b}\|$', 'FontSize', 20, 'Interpreter', 'latex')
legend('GMRES', 'IDR(8)', 'BiCGSTAB', 'FontSize', 15, 'Interpreter', 'latex')