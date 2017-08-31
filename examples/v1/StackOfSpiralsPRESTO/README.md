# TOPPE sequence example 

## 3D stack-of-spirals PRESTO functional MRI sequence

Main script is 'main.m'. See ../README.md for further explanation.

### Preview sequence in MATLAB
```
>> addpath('../..');
>> playseq('scanloop.txt',2,1,0.03);
```

### Acquisitions parameters

+ dynamic 3.33mm iso
+ Volume TR is about 3sec
+ 4 fully sampled frames at beginning, then a series of R=6 undersampled (3 in-plane x 2 in kz) fMRI frames.
