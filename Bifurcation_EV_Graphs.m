%Eigenvalue and bifurcation graphs

function Bifurcation_EV_Graphs()
gamma = 1;
%alpha = 0.29
%alpha = 2.347510026817774; %transition time from low to high density; not infectivity!
beta = 3.374; %calculated from root-finding the imaginary values in the Jacobian
mu = 0.01;
Lambda = 1;

%R0 = (Lambda/mu)*(beta/gamma); %here we assume a steady population where deaths and immigration over a three-week interval are negligible
%R0 = math.sqrt(Lambda/mu)*math.sqrt(beta/gamma) #here we assume a steady population where deaths and immigration over a three-week interval are negligible
%printf(R0)

alpha_min = 0.001;
alpha_max = 5;

alpha_evaluate = linspace(alpha_min, alpha_max, 500);
%transient_time = 14; %model peaks quickly

bifurcation_alpha = zeros(length(alpha_evaluate));
bifurcation_I = zeros(length(alpha_evaluate));
bifurcation_S = zeros(length(alpha_evaluate));
bifurcation_E = zeros(length(alpha_evaluate));
bifurcation_R = zeros(length(alpha_evaluate));

ev4 = zeros(length(alpha_evaluate));
ev1 = zeros(length(alpha_evaluate));
ev2 = zeros(length(alpha_evaluate));
ev3 = zeros(length(alpha_evaluate));




for val = (1:length(alpha_evaluate))
    alpha = alpha_evaluate(val);
    I_steady = ((Lambda*alpha)/((gamma+mu)*(alpha+mu)))-(mu/beta);
    S_steady = ((gamma+mu)*(alpha+mu))/(alpha*beta);
    E_steady = gamma/(alpha+mu) - (mu*(gamma+mu))/(alpha*beta);
    R_steady = (gamma*Lambda*alpha)/(mu*(gamma+mu)*(alpha+mu)) - gamma/beta;
    bifurcation_I(val) = I_steady;
    bifurcation_S(val) = S_steady;
    bifurcation_E(val) = E_steady;
    bifurcation_R(val) = R_steady;

    bifurcation_alpha(val) = alpha;

    Jac = [-(Lambda*alpha*beta)/((alpha+mu)*(gamma+mu))-mu, 0, -((gamma+mu)*(alpha+mu))/alpha, 0; 
    (Lambda*alpha*beta)/((alpha+mu)*(gamma+mu)), -alpha-mu, ((gamma+mu)*(alpha+mu))/alpha, 0;
    0, alpha, -gamma-mu, 0;
    0, 0, gamma, mu];

    eigs = eig(Jac);
    ev1(val) = eigs(1);
    ev2(val) = eigs(2);
    ev3(val) = eigs(3);
    ev4(val) = eigs(4);
end    

bif = [bifurcation_S; bifurcation_I; bifurcation_E; bifurcation_R]';

figure;
    plot(bifurcation_alpha, bifurcation_S, 'LineWidth', 2);
    hold on;
    plot(bifurcation_alpha, bifurcation_E, 'LineWidth', 2);
    plot(bifurcation_alpha, bifurcation_I, 'LineWidth', 2);
    plot(bifurcation_alpha, bifurcation_R, 'LineWidth', 2);
    hold off;
    xlabel('Progression Rate (\alpha)');
    ylabel('Number per Compartment');
    set(gca, 'YScale', 'log');
    legend('Susceptible (S)', 'Exposed (E)', 'Infectious (I)', 'Recovered (R)');
    title('Steady State Bifurcation in MATLAB');
    grid on;

    im_ev1 = real(ev1);
    im_ev2 = real(ev2);
    im_ev3 = real(ev3);
    im_ev4 = real(ev4);

    figure;
    plot(bifurcation_alpha, im_ev1, 'LineWidth', 2);
    hold on;
    plot(bifurcation_alpha, im_ev2, 'LineWidth', 2);
    plot(bifurcation_alpha, im_ev3, 'LineWidth', 2);
    plot(bifurcation_alpha, im_ev4, 'LineWidth', 2);
    hold off;
    xlabel('Progression Rate (\alpha)');
    ylabel('Real Part of Eigenvalues');
    %set(gca, 'YScale', 'log');
    legend('1', '2', '3', '4');
    title('Real Eigenvalues of Jacobian in MATLAB');
    grid on;

im_ev1 = imag(ev1);
    im_ev2 = imag(ev2);
    im_ev3 = imag(ev3);
    im_ev4 = imag(ev4);

    figure;
    plot(bifurcation_alpha, im_ev1, 'LineWidth', 2);
    hold on;
    plot(bifurcation_alpha, im_ev2, 'LineWidth', 2);
    plot(bifurcation_alpha, im_ev3, 'LineWidth', 2);
    plot(bifurcation_alpha, im_ev4, 'LineWidth', 2);
    hold off;
    xlabel('Progression Rate (\alpha)');
    ylabel('Imaginary Part of Eigenvalues');
    %set(gca, 'YScale', 'log');
    legend('1', '2', '3', '4');
    title('Imaginary Eigenvalues of Jacobian in MATLAB');
    grid on;


end