% Tests loadpfile vs loadpfilefast
clear; clc; close all;
pfile = 'Ptest1.7';
%pfile = 'Ptest2.7';

%% Test old file
tic
dat1 = toppe.utils.loadpfile(pfile);
t1 = toc;


%% New one
tic
dat2 = toppe.utils.loadpfilefast(pfile);
t2 = toc;

% Check that data is the same
if all(dat1(:)==dat2(:))
    disp('Outputs match!');
else
    error('Output mismatch.');
end

% Compare times
fprintf('loadpfile time: %0.1f seconds\n',t1);
fprintf('loadpfilefast time: %0.1f seconds\n',t2);
fprintf('Speedup factor: %0.3f\n',t1/t2);