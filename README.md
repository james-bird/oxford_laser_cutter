# oxford_laser_cutter

Code for generating cutting paths for A Series Oxford Laser micro machining system from .dxf files. May work for other Oxford Laser systems too, but this is untested.

Addresses the need for cutting large items easily using the high cut speeds from the galvos. 

The code takes a .dfx file and breaks it up into overlapping squares which are then moved to with the bed, and then cut with the galvos. outputs a .pgm program file to acheive this.

To run: Set cutting and file paramters in dxf2lincut.m, and then run the script. Plots of cut-paths are produced to aid any debugging. Additionally, a scaling variable has been added so cutting paths can be massively reduced in size for rapid testing.

The script is (very) unoptimised, but I no-longer have access to the machine to reliably make changes. 

The reading of the .dxf was written for .dxf files produced with CREO. If they are generated with other software, it is likely you might have to make some small changes to the parsing in the code.

This script saved me countless hours, and so might be of use to other people working with Oxford Lasers micromachine systems!

Use at your own risk. 
