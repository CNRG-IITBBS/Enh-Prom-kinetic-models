clear all;

%%% Parameters %%%%%%%%
f = 1.0;
RC_values = 5.0;
enhancers=20;
for RC=RC_values
    for num_enhancers = enhancers
        num_molecules = num_enhancers + 2;
        slopes_1_values = 0.99;
        slopes_2_values = 1.01;
       
        
        for s1 = slopes_1_values
            for s2 = slopes_2_values
                slopes(1) = s1;
                slopes(2) = s2;
                r_values = 0.4;
                g_values = 0.4;
                for r = r_values
                    for g = g_values
                        ma_values =5.0;
                        for ma = ma_values
                            mr_values = 5.0;
                            for mr = mr_values
                                for s = 1:20
                                    runs = 1000000;
    
                                    %%% Defining state %%%%
                                    enh = zeros(3, num_molecules);
                                    rate = cell(3, num_molecules);
    
                                    for i = 1:size(rate, 1)
                                        for j = 2:size(rate, 2) - 1
                                            rate{1, j} = zeros(2); % Initialize rate matrix for activator
                                            rate{2, j} = zeros(2); % Initialize rate matrix for enhancer
                                            rate{3, j} = zeros(2); % Initialize rate matrix for repressor
                                        end
                                    end
    
                                    enh(2, 2) = 1; % Initial state
                                    time = 0;
                                   % fileID1 = fopen(sprintf('enhancers_matrix_%d_r_%.5f_g_%.5f_s1_%.3f_s2_%.3f_ma_%.2f_mr_%.2f_RC_%.2f_AC5_t%d.txt',num_enhancers, r, g, slopes(1), slopes(2), ma, mr,RC, s), 'w');
                                    fileID2 = fopen(sprintf('enh_num_%d_r_%.5f_g_%.5f_s1_%.3f_s2_%.3f_ma_%.2f_mr_%.2f_RC_%.2f_AC5_t%d.txt', num_enhancers, r, g, slopes(1), slopes(2), ma, mr,RC, s), 'w');
                                    %fileID2=fopen('cluster.txt','w');
                                    %fileID2 = fopen(sprintf('enh_num_%d_s1_%d_s2_%d_ma_%d_mr_%d_AC_100_t%d.txt', num_enhancers,slopes(1),slopes(2),ma,mr,s), 'w');
                                    avg_col = 0.0;
                                    
                                    for i = 1:runs
                                        % Generate new rate values
                                        AC = rect(i); % This could be updated based on your specific logic
                                        [a, b, c, d, u, v, e_to_ea, ea_to_e, e_to_er, er_to_e] = generate_values(slopes, num_molecules, r, g, ma, mr, AC, RC, f);
                                        %disp(a);
                                        %disp(b);
                                        %disp(c);
                                        %disp(d);
                                        %disp(u);
                                        %disp(v);
                                        %disp(e_to_er);
                                        %disp(er_to_e);
                                        %disp(e_to_ea);
                                        %disp(ea_to_e);
                                        % Update rate matrices
                                        for j = 2:size(rate, 2) - 1
                                            rate{1, j}(1, 1) = 0; % up rate
                                            rate{1, j}(1, 2) = ea_to_e(j); % down rate
                                            rate{1, j}(2, 1) = c(j); % right rate
                                            rate{1, j}(2, 2) = d(j - 1); % left rate
    
                                            rate{2, j}(1, 1) = e_to_ea(j); % up rate
                                            rate{2, j}(1, 2) = e_to_er(j); % down rate
                                            rate{2, j}(2, 1) = a(j); % right rate
                                            rate{2, j}(2, 2) = b(j - 1); % left rate
    
                                            rate{3, j}(1, 1) = er_to_e(j); % up rate
                                            rate{3, j}(1, 2) = 0; % down rate
                                            rate{3, j}(2, 1) = u(j); % right rate
                                            rate{3, j}(2, 2) = v(j - 1); % left rate
                                        end
    
                                        % Simulation
                                        [row, col] = find(enh == 1);
                                        index_i = row;
                                        index_j = col;
    
                                        %avg_col = avg_col + col*t;
    
                                        c1 = rate{index_i, index_j}(1, 1); % up
                                        c2 = rate{index_i, index_j}(1, 2); % down
                                        c3 = rate{index_i, index_j}(2, 1); % right
                                        c4 = rate{index_i, index_j}(2, 2); % left
    
                                        R = c1 + c2 + c3 + c4; % sum rates of all the reactions
                                        r2 = R * rand(); % generate random number between 0 and R
                                        r1 = rand(); % another random number for time
                                        t = (1 / R) * log(1 / r1); % time from exponential distribution
                                        time = time + t; % time at which next event happens
                                        avg_col = avg_col + (col-1)*t;
                                        % Update enhancer states based on reaction rates
                                        if index_i > 0 && r2 > 0 && r2 <= c1
                                            enh(index_i, index_j) = 0;
                                            enh(index_i - 1, index_j) = 1;
                                        elseif index_i < 3 && r2 > c1 && r2 <= c1 + c2
                                            enh(index_i, index_j) = 0;
                                            enh(index_i + 1, index_j) = 1;
                                        elseif index_j < num_molecules && r2 > c1 + c2 && r2 <= c1 + c2 + c3
                                            enh(index_i, index_j) = 0;
                                            enh(index_i, index_j + 1) = 1;
                                        elseif index_j > 1 && r2 > c1 + c2 + c3 && r2 <= c1 + c2 + c3 + c4
                                            enh(index_i, index_j) = 0;
                                            enh(index_i, index_j - 1) = 1;
                                        end
                                        
                                        % Logging the state
                                        %   fprintf(fileID1, 'Step %f:\n', time);
                                        %  for row = 1:size(enh, 1)
                                        %      fprintf(fileID1, '%f ', enh(row, :));
                                        %     fprintf(fileID1, '\n');
                                        %  end
                                        % fprintf(fileID1, '\n');
    
    
                                        fprintf(fileID2, '%f %f %f %d\n', time, AC, t, col-1);
                                         %disp(i);disp(col);disp(t);disp(time);
                                       %disp(enh);
                                        %disp(avg_col);
                                    end
                                    
                                    cl = avg_col / time;
                                    %disp(cl);
    
    %                                 fclose(fileID1);
                                    fclose(fileID2);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
%%% Function to generate the rates
function [a, b, c, d, u, v, e_to_ea, ea_to_e, e_to_er, er_to_e] = generate_values(slopes, num_molecules, r, g, ma, mr, AC, RC, f)
    % Initialize arrays for each variable with initial values
    a = zeros(1, num_molecules);
    b = zeros(1, num_molecules);
    c = zeros(1, num_molecules);
    d = zeros(1, num_molecules);
    u = zeros(1, num_molecules);
    v = zeros(1, num_molecules);
    e_to_ea = zeros(1, num_molecules);
    ea_to_e = zeros(1, num_molecules);
    e_to_er = zeros(1, num_molecules);
    er_to_e = zeros(1, num_molecules);
    for i = 2:2
        a(i) = f* r;
        b(i) = f* g;
    end
    for i = 3:(num_molecules - 2)
        a(i) = f* max(0, r * ((slopes(1)^(i - 2))));
        b(i) = f* max(0, g * ( (slopes(2)^(i - 2) )));
    end
    for i = 2:2
        c(i) = f * ma * r;
        d(i) = f * g;
        u(i) = f*r;
        v(i) = f* mr * g;
    end
    for i = 3:(num_molecules - 2)
        c(i) = f * max(0, (ma * r) * ((slopes(1)^(i - 2) )));
        d(i) = f * max(0, g * ( slopes(2)^ (i - 2)));
        u(i) = f* max(0, (r) * ((slopes(1)^(i - 2))));
        v(i) = f* max(0, (mr * g) * ((slopes(2)^(i - 2))));
    end

    for i = 2:(num_molecules - 1)
        e_to_ea(i) = f * AC * r;
        ea_to_e(i) = f * g;
        e_to_er(i) = f* RC * r;
        er_to_e(i) = f * g;
    end
end
%function AC = rect(i)
%    T = 100000; % Cycle duration
%    amplitude = 1.1; % Amplitude of the rectangular pulse
%    AC = amplitude * (mod(i, T) > T / 2);
%end

function AC = rect(i)
    T = 100000;  % Total cycle duration
    on = 5.0;  
    off = 0.0; 

    % Each "on" and "off" state has a width of 0.1 * T
    width1=0.5;
    width = width1 * T;
    num_stages = 1/width1;  % Total number of stages (5 "on" and 5 "off")

    % Calculate the total duration of the 10 alternating states
    total_cycle_duration = num_stages * width;

    % Determine the position within the total cycle
    position_in_cycle = mod(i, total_cycle_duration);

    % Determine which stage (on/off) the current position falls into
    stage = floor(position_in_cycle / width) + 1;

    % Alternate between "on" and "off" states
    if mod(stage, 2) == 0
        AC = on;  % "On" state for odd stages
    else
        AC = off; % "Off" state for even stages
    end
end
