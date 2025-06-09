EMDLAB finite element package

Electrical Machines Design Laboratory

Overview:
EMDLAB is an open-source numerical package developed in the MATLAB environment for the design
and analysis of electrical machines, including motors, generators, transformers, and actuators.

The package is organized as a library of MATLAB objects. Depending on your specific design needs, 
you can select the appropriate modules to obtain your desired results. With EMDLAB, 
you can also develop customized standalone software tailored to your applications.

Setup Instructions (Windows 64-bit):
1) Download the emdlab-win64.zip file.
2) Extract the zip file and place "emdlab-win64" folder to the "C:\" directory without changing the folder name.
3) To use EMDLAB package in your MATLAB code, add following line at the begining of your mfile:
--->> addpath(genpath('C:\emdlab-win64'));
