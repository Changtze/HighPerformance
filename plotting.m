%%% test case 6 %%%

case1 = 'D:\Aeronautics\HPC\coursework\case1\energy.txt';

case2 = "D:\Aeronautics\HPC\coursework\case2\energy.txt";

case3_e = "D:\Aeronautics\HPC\coursework\case3\energy.txt";
case3_o = "D:\Aeronautics\HPC\coursework\case3\output.txt";

case4_e = "D:\Aeronautics\HPC\coursework\case4\energy.txt";
case4_o = "D:\Aeronautics\HPC\coursework\case4\output.txt";

case5_e = "D:\Aeronautics\HPC\coursework\case5\energy.txt";
case5_o = "D:\Aeronautics\HPC\coursework\case5\output.txt";

case6_e = "D:\Aeronautics\HPC\coursework\case6\energy.txt";
case6_o = 'D:\Aeronautics\HPC\coursework\case6\output.txt';

% test case 1 
opts_1 = detectImportOptions(case1, 'Delimiter', ' ', 'VariableNamesLine', 1);
energy_1 = readtable(case1, opts_1);


% test case 2 
opts_2 = detectImportOptions(case2, 'Delimiter', ' ', 'VariableNamesLine', 1);
energy_2 = readtable(case2, opts_2);

% test case 3 
opts_3 = detectImportOptions(case3_o, 'Delimiter', ' ', 'VariableNamesLine', 1);
opts_3_e = detectImportOptions(case3_e, 'Delimiter', '  ', 'VariableNamesLine', 1);

output_3 = readtable(case3_o, opts_3);
energy_3 = readtable(case3_e, opts_3_e);

% test case 4 
opts_4 = detectImportOptions(case4_o, 'Delimiter', ' ', 'VariableNamesLine', 1);
opts_4_e = detectImportOptions(case4_e, 'Delimiter', ' ', 'VariableNamesLine', 1);
output_4 = readtable(case4_o, opts_4);
energy_4 = readtable(case4_e, opts_4_e);

% test case 5 
opts_5 = detectImportOptions(case5_o, 'Delimiter', ' ', 'VariableNamesLine', 1);
opts_5_e = detectImportOptions(case5_e, 'Delimiter', ' ', 'VariableNamesLine', 1);
output_5 = readtable(case5_o, opts_5);
energy_5 = readtable(case5_e, opts_5_e);

% test case 6
opts_6 = detectImportOptions(case6_o, 'Delimiter', ' ', 'VariableNamesLine', 1);
opts_6_e = detectImportOptions(case6_e, 'Delimiter', ' ', 'VariableNamesLine', 1);
output_6 = readtable(case6_o, opts_6);
energy_6 = readtable(case6_e, opts_6_e);


% % Plot
% plot(energy_1.t, energy_1.E)
% xlabel("Time")
% ylabel("Kinetic energy")
% title("Test case 1: Kinetic Energy")
% hold off
% 
% plot(energy_2.t, energy_2.E)
% xlabel("Time")
% ylabel("Kinetic energy")
% title("Test case 2: Kinetic Energy")
% 
% 
% plot(energy_3.t, energy_3.E)
% xlabel("Time")
% ylabel("Kinetic energy")
% title("Test case 3: Kinetic Energy")

% min_particle_distance(output_3)
% plot(output_3.x(output_3.pNum == 1), output_3.y(output_3.pNum == 1))
% hold on
plot_particle_distance(output_3, 1, 2)
% plot(output_3.x(output_3.pNum == 2), output_3.y(output_3.pNum == 2))
% xlabel("x")
% ylabel("y")
% title("Test case 3: Trajectory")

