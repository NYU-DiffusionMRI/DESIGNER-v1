function processStudy(studyDir);
%   processStudy is a MATLAB function designed to provide "hands-off"
%   preprocesxing to diffusion MRI images.
%
%   processStudy(studyDir) procprocesses the input directory 'studyDir' and
%   returns diffusion/kurtosis maps
%
%   Author: Siddhartha Dhiman
%   Date Created: 02/14/2019
%   Matlab 2018b

%% Start Fresh
%   Clear everything to minimize errors
clc; close all; clear all; warning off;

%% Locate Study Directory
%   Opens GUI to locate directory
studyPath = uigetdir(pwd,'Select subject directory');
studyDir = dir(studyPath);
studyDirFolders = dir(studyPath);

%%
