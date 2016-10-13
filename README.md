# analysis
Code to analyze multipopout, priority, and delayed task data

**readall.m**
* read .plx files and convert usable data to matlab variables stored in a .mat files
* need to run this before any analysis code will work. I tried to incorportate it into my analysis code so that if the proper .mat file isn't found, it'll automatically run readall.m

**runles_working.m**
* read file names from excel and analyze delayed task data
* outputs to same excel file if there's a visual burst and a motor burst

**listTrials.m**
* read file names from excel and formats eyedata and other info to a .mat for use with Manual_MS_correction.m

**Manual_MS_correction.m**
* open a GUI to manually decide if a correct saccade happened

**batch.m**
* run series of analysis on main tasks
