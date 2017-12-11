# Pulseq to TOPPE conversion

Disclaimer: this conversion may require manual tweaking.

From a .seq file, generate the following files:

	.mod files:     one for each module

	modules.txt:    list of .mod files

	scanloop.txt:   defines pulse sequence (sequence of modules; waveform amplitudes; RF/ACQ phase; etc)


### Run the example

In Matlab, do
```
>> addpath ~/github/toppe/matlab/lib/v1    % edit accordingly
>> main;
```



