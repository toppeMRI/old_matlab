# TOPPE sequence example 

## 3D RF-spoiled gradient-recalled echo (SPGR/FLASH/T1-FFE)

This folder contains the MATLAB script **'main.m'** which generates the file **'scan.tgz'** containing the following TOPPE files:

+ modules.txt
+ scanloop.txt
+ tipdown.mod
+ readout.mod
+ spoiler.mod

In MATLAB, simply do:

``` >> main; ```

To **scan**, copy scan.tgz to /usr/g/bin/ on the scanner, open a console and type:

```
$ cd /usr/g/bin/;
$ tar xzf scan.tgz
```

Next, prescribe and run the toppe sequence (currently named 'toppe9a').


