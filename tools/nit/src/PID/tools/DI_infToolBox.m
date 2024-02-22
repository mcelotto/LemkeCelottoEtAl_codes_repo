function [ di, diSh ] = DI_infToolBox( X, Y, hY, bias, computeShuffled )

% Computing DI using infToolBox
% DI(X,Y) = I(hY,X ; Y) - I(hY, Y);

% I(hY,X ; Y)
%[R, nt] = buildr(Y, [hY; X]);

diSh = -1;

%Opt.nt = nt;
Opt.method = 'dr';
Opt.bias   = bias;
Opt.trperm = 0;
Opt.bin_method = 'none';
Opt.n_bins = 2;

if computeShuffled
    [IhyXY, IhyXYsh] = information(R, Opt, 'I', 'Ish');
    
    % I(hY; Y)
    [R, nt] = buildr(Y, hY);
    Opt.nt = nt;

    [IhyY, IhyYsh] = information(R, Opt, 'I', 'Ish');

    di = IhyXY - IhyY;
    diSh = IhyXYsh - IhyYsh;
else
    IhyXY = information([hY; X], Y, Opt, {'I'});

    % I(hY; Y)
    %[R, nt] = buildr(Y, hY);
    %Opt.nt = nt;

    IhyY = information(hY,Y, Opt, {'I'});

    di = IhyXY{1} - IhyY{1};
end
end

