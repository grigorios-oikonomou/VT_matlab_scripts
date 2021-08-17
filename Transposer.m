% Column breaker
powers = 14 % the # of powers that you are using
trials = 30 % the # of trials for each power (has to be the same for all powers)

for i=1:powers
    b(:,i) = a(((i-1)*trials+1):((i-1)*trials+trials))
end

% you need to check the above for outliers at prism, and what you get as
% clean run into the transposer


%% Transposer
powers = 14 % the # of powers that you are using
trials = 30 % the # of trials for each power (has to be the same for all powers)
% a is the column of MovedRatio
% b is the transpose of the MovedRatio that you can paste into Prism

for i=1:powers
    b(i,:)=a(((i-1)*trials+1):((i-1)*trials+trials))'
end


%% matrix into column from smaller to larger
b = sort(a(:))

%% array into single row with column after column
b = a(:)';

%% array into single column with column after column
b = a(:)


%% array into single column with row after row
b = a'
c = b(:)


sort(c)
%% positions of 0 and non-0

d = find(c==0)
e = find(c>0)

%% all above 0 =1
