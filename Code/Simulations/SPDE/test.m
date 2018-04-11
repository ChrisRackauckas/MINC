id = timeStamp();
K = 10000;
parforTimer = tic;
parfor_progressf(id,parforTimer,K);
parfor i=1:K
    pause(10);
    parfor_progressf(id,parforTimer);
end
parfor_progressf(id,parforTimer,0);
toc(parforTimer);