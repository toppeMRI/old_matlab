
ge2seq.m: 
Convert TOPPE files to a Pulseq (.seq) file.

Uses Pulseq Matlab library that can be downloaded as follows:
```
 $ git clone https://github.com/pulseq/pulseq/
```

Important: uncomment the line in makeArbitraryRf (Pulseq) that scales the b1, since
rf2pulseq() is already scaled correctly (in units of Hz).

Usage examples:
```
>> ge2seq();    % assumes scan files are 'scanloop.txt', 'modules.txt', and 'timing.txt'
>> ge2seq('timingfile', 'newtiming.txt');
```
