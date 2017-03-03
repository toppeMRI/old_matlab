# TOPPE sequence examples

Working with TOPPE involves three basic steps:

1. Create the required sequence files: modules.txt, scanloop.txt, and the .mod files.
2. Pre-view the sequence using playseq.m.
3. Scan.


## Step 1: Create TOPPE sequence files

Each subfolder contains a MATLAB script **'main.m'** which generates the file **'scan.tgz'**, a tar file that contains the following TOPPE files:

+ modules.txt
  + Lists the .mod files, and indicates whether a module is an RF excitation, readout, or gradients-only module.
+ scanloop.txt
  + Lists the order in which to play out the modules, as well as all other dynamic sequence information.
+ tipdown.mod (a .mod file with this name *must* exist on the scanner)
  + RF excitation module that excites a 2cm slab.
+ readout.mod (a .mod file with this name *must* exist on the scanner)
  + 3D Cartesian (spin-warp) readout module.
+ Any additional .mod files listed in modules.txt.

In MATLAB, simply do:

``` >> main; ```


## Step 2: Pre-view the sequence in MATLAB

Use the script playseq.m to view the sequence in movie loop mode, e.g.,:

```
>> playseq('scanloop.txt',nModPerTR,nTRskip);
```
where nModPerTR is the number of modules per TR, and nTRskip is the number of TRs to skip (for faster looping).


## Step 3: Scan

To **scan**, copy scan.tgz to /usr/g/bin/ on the scanner, open a console and type:

```
$ cd /usr/g/bin/;
$ tar xzf scan.tgz
```

Next, prescribe and run the toppe sequence (e.g., 'toppev1').



