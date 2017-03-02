# TOPPE sequence example 

## 3D RF-spoiled gradient-recalled echo (SPGR/FLASH/T1-FFE)

In MATLAB, do:
>> main;

This will generate scanloop.txt, and create scan.tgz containing TOPPE sequence files for this scan.

To scan, copy scan.tgz to /usr/g/bin/ on the scanner, open a console and type:
>> cd /usr/g/bin/;
>> tar xzf scan.tgz
Then, prescribe and run the toppe9a sequence.


