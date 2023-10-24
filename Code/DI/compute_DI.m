function [info] = compute_DI(hX, hY, Y, opts)

%%% Script to compute Directed Information (DI) from the emitter X to the receiver Y

%%% - *hX*: must be an array of *nDimsX X nTrials* response matrix describing the past of the emitter variables on each of the *nDims* dimensions for each trial.
%%% - *hY*: must be an array of *nDimsY X nTrials* response matrix describing the past of the receiver variable on each of the *nDims* dimensions for each trial.
%%% - *Y*: must be an array of *nDimsY X nTrials* response matrix describing the response of each of the *nDims* dimensions for each trial.
%%% - *opts*: options used to calculate II (see further notes).

I1 = information(hY,Y,opts,{'I'});
I2 = information([hY; hX],Y,opts,{'I'});

info = I2{1}(1)-I1{1}(1);

end