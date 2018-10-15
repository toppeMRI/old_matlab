# Matlab scripts for TOPPE v2

These scripts work with the toppev2.e driver/psd.

Changes/fixes from toppev1:

 * Fixed theta waveform playback bug (especially important for, e.g., tailored RF pulses or high time-bandwidth SLR pulses)
 * It's now possible to associate multiple waveforms with each .mod file; accordingly, an additional column is required in scanloop.txt specifying which of the waveform(s) to play out.
 * The header in the .mod file has changed, so be sure to use the mat2mod.m in this folder when scanning with toppev2. To convert a .mod file from v1 to v2, you may use the script v1tov2.m in this folder.
 * Changed the interface/usage of most of the .m functions; now uses the ('param', value) syntax.
 * mat2mod.m now includes checks for system hardware limits, and should be more robust.
 * getScanTime.m: returns total scan time (new)


 



