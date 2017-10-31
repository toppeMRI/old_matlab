# Basic TOPPE MATLAB functions

Matlab scripts for TOPPE v2.

Changes/fixes from toppev1:

 * fixed theta waveform playback bug (especially important for, e.g., tailored RF pulses or high time-bandwidth SLR pulses)
 * it's now possible to associate multiple waveforms with each .mod file; accordingly, an additional column is required in scanloop.txt specifying which of the waveform(s) to play out
 * A TOPPE to Pulseq conversion script, ge2seq.m, is now included in the Matlab library
 * The header in the .mod file has changed, so be sure to use the mat2mod.m in this folder when scanning with toppev2

