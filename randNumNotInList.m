function [rNum] = randNumNotInList_old(N, usedList)



A = 1:N;

A(ismember(A,usedList)) = [];

msize = numel(A);

rNum = A(randperm(msize,1));


end