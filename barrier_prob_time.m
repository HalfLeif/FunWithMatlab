%% Methods for plotting

ftitle = @(s) title(s, 'FontSize', 16);
fxlabel = @(s) xlabel(s, 'FontSize', 14);
fylabel = @(s) ylabel(s, 'FontSize', 14);
setlegend = @(hl) set(hl, 'FontSize', 14);
flegend = @(c) setlegend(legend(c{:}));


%% Bernoulli version
f = @(n, p) 1 - (1 - p).^n;
len = @(x) size(x, 2);

% p = prob of long turn
step = 0.20;
p = (1 - step):-step:0.01;

maxn = 7;
processes = 0:0.2:maxn;
num_lines = len(processes);

out = zeros(num_lines, size(p, 2));
for i = 1:num_lines
    out(i, :) = f(processes(i), p);
end

plot(processes, out)
ftitle('Bernoulli model')
fxlabel('num processes/players')
fylabel('estimated time')
legend_cell = cellstr(num2str(p', 'p=%.2f'));
flegend(legend_cell)

clear step maxn f


%% Exp version

e_cdf = @(x, lambda) 1 - exp(-lambda .* x);
e_pdf = @(x, lambda) lambda .* exp(-lambda .* x);

% Supports vector argument for x
total_pdf = @(n, x, lambda) ...
            n .* e_pdf(x, lambda) .* (e_cdf(x, lambda) .^ (n - 1));

for i = 1:num_lines
    n = processes(i);
    for j = 1:len(p)
        lambda = 1 ./ p(j);
        integ = @(x) x .* total_pdf(n, x, lambda);
        out(i, j) = integral(integ, 0, 1000 * lambda);
    end
end

figure
plot(processes, out)
ftitle('Exp model, \lambda = 1/p')
fxlabel('num processes/players')
fylabel('estimated time')
legend_cell = cellstr(num2str(p', 'p=%.2f'));
flegend(legend_cell)

clear lambda integ out j


%% Plot pdf for different num players

figure
hold on

tp = 0.5;
players = 1:2:7;
xs = 0:0.02:4;
legend_cell = cell(1, len(players));
for i = 1:len(players)
    n = players(i);
    pdf = @(x) total_pdf(n, x, 1 ./ tp);
    plot(xs, pdf(xs))
    legend_cell(i) = {[num2str(n,'players=%i'),' ', num2str(tp,'p=%.2f')]};
end
ftitle('Pdf for different number of players/processes')
flegend(legend_cell)
fxlabel('time')
fylabel('probability')
hold off
clear xs pls tp i n pdf legendcell