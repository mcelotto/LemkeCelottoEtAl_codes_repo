disp('Testing GLMnet:')

x=randn(100,20);
g4=randsample(4,100,true);
fit3=glmnet(x,g4,'multinomial');
opts=struct('mtype','grouped');
fit3a=glmnet(x,g4,'multinomial',opts);

disp('Test OK');