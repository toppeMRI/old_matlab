# TOPPE sequence example 

## 3D RF-spoiled gradient-recalled echo (SPGR/FLASH/T1-FFE)

This folder contains the MATLAB script **'main.m'** which generates the file **'scan.tgz'** containing the following TOPPE files:

+ modules.txt
  + lists the three .mod files, and indicates whether a module is an RF excitation, readout, or gradients-only module.
+ scanloop.txt
  + Lists the order in which to play out the modules, as well as all other dynamic sequence information.
+ tipdown.mod (a .mod file with this name *must* exist on the scanner)
  + RF excitation module that excites a 2cm slab.
+ readout.mod (a .mod file with this name *must* exist on the scanner)
  + 3D Cartesian (spin-warp) readout module.
+ spoiler.mod
  + A z gradient trapezoid.

In MATLAB, simply do:

``` >> main; ```

To **scan**, copy scan.tgz to /usr/g/bin/ on the scanner, open a console and type:

```
$ cd /usr/g/bin/;
$ tar xzf scan.tgz
```

Next, prescribe and run the toppe sequence (currently named 'toppe9a').

### Acquisitions parameters

+ dynamic 3.33mm iso
+ TR is about 18e-3*(3*60+nskip) = ~3sec
+ R=6 undersampling (3 in-plane x 2 in kz)
