%% Erik Rutyna
close all
clear
clc

N = 32;
% Re = 100;
% Re = 400;
Re = 1000;
h = 1 / N;

% These two functions were written by Jose Luiz Vargas de
% Mendonca. All they do is return the  list of velocities from
% the paper for a given Re. He was kind enough to share it so
% that we don't make typos when translating the numbers from text
% to code in MATLAB. I hope using them is ok.
[yGhia, uGhai] = U_vel(Re);
[xGhia, vGhia] = V_vel(Re);

rFile = ['residuals' '.csv'];
uFile = ['uVelocity' '.csv'];
vFile = ['vVelocity' '.csv'];

% rFile = ['errorNorm_' num2str(N) '_' num2str(Re) '.csv'];
% uFile = ['uVelocity_' num2str(N)  '_' num2str(Re) '.csv'];
% vFile = ['vVelocity_' num2str(N) '_' num2str(Re) '.csv'];

residuals = readtable(rFile);
uVel = readtable(uFile);
vVel = readtable(vFile);

residuals = table2array(residuals);
uVel = table2array(uVel(:, 1:N+3));
vVel = table2array(vVel(:, 1:N+2));

uVel = uVel(2:33, 2:34);
vVel = vVel(2:34, 2:33);

titleS = ['N = ' num2str(N) ' , Re = ' num2str(Re), ];
uYlabel = ['x-velocity'];
vYlabel = ['y-velocity'];

figure(1)
hold on
scatter(yGhia, uGhai)
plot(linspace(h/2, 1-h/2, N), uVel(:, 17))
xlabel('y-Coordinate', 'interpreter', 'latex')
ylabel(uYlabel, 'interpreter', 'latex')
% title(titleS, 'interpreter', 'latex')
legend({'Ghia-Shin','My simulation'}, 'location', 'northwest')
grid on

figure(2)
hold on
scatter(xGhia, vGhia)
plot(linspace(h/2, 1-h/2, N), vVel(17, :))
xlabel('x-Coordinate', 'interpreter', 'latex')
ylabel(vYlabel, 'interpreter', 'latex')
legend({'Ghia-Shin','My simulation'})
% title(titleS, 'interpreter', 'latex')
grid on


levels = [-0.1175, -0.1150, -0.1100, -0.100, -0.09, -0.07, -0.05 -0.03, -0.01, -1e-4, -1e-5, -1e-7, -1e-10...
                1e-8 1e-7 1e-6 1e-5 5e-5 1e-4 2.5e-4 5e-4 1e-3 1.5e-3 3e-3];
            
psi = zeros(N+1, N);

for i=2:N            
    psi(:, i) = psi(:, i-1) - vVel(:, i) / N;
end

X = linspace(h/2, 1-h/2, N);
Y = linspace(0, 1, N+1);
[X, Y] = meshgrid(X, Y);

figure(3)
contourf(X, Y, psi, levels)
colormap(jet)
a = colorbar;
a.Label.String = 'Stream Function (\Psi)';
xlabel('x-Coordinate', 'interpreter', 'latex')
ylabel('y-Coordinate', 'interpreter', 'latex')
% title(titleS, 'interpreter', 'latex')


figure(4)
plot([1:length(residuals)], residuals)
xlabel('Iteration number', 'interpreter', 'latex')
ylabel('Residual', 'interpreter', 'latex')
set(gca, 'YScale', 'log')
% title(titleS, 'interpreter', 'latex')
grid on