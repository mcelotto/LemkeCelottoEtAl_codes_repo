%check if infotoolbox is installed correctly
clear all
clc
close all

disp("Testing MI - Binned Methods:")

S = ones(1,100);
S(51:100) = 2;
R = ones(1,100);
R(1:10) = 1:10;
R(11:30) = 11:30;
R(31:100) = 30;
opts.verbose = true;
methods = ["dr", "gs"];
biases = ["qe", "pt", "gsb", "naive"];

disp("Estimation methods:")
disp("--------------------------");
disp("| 'dr' | Direct method   |");
disp("| 'gs' | Gaussian method |");
disp("--------------------------");
disp("Bias reduction methods:")
disp("-------------------------------------");
disp("| 'qe'    | Quadratic EXtrapolation |");
disp("| 'pt'    | Panzeri & Treves 1996   |");
disp("| 'gsb'   | Gaussian bias           |");
disp("| 'naive' | Biased naive estimates  |");
disp("-------------------------------------");

for m = 1:length(methods)
    for b = 1:length(biases)
        bias = biases(b);
        method = methods(m);
        if (strcmp(bias,'gsb') && ~strcmp(method,'gs')) || (strcmp(bias,'pt') && ~strcmp(method,'dr'))
            % bias method gsb can only be applied with gs method
            % bias method pt can only applied with dr method so all
            % other combinations are skipped
            continue
        end
        disp(['Testing method ', method, ' with bias reduction ', bias]);
        opts.method = method;
        opts.bias = bias;
        opts.btsp = 100;
        opts.bin_methodX = 'eqpop';
        opts.n_binsX = 3;
        outputs = information(R, S, opts, {'I'});
    end
end

disp('Test OK');
