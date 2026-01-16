% Convert the cell array to a table for better handling and export
clear all;close all;clc

allUCBCasesResults = table();
sleepStatus ='Sleep';
phaseOfStudy = 'Pre';
subjectStatus = 'case';


load('Sleep1_Pre_results.mat')
[allUCBCasesResults ] = get_cvs_from_computational_metrics_mat(allUCBCasesResults,eegComputationalMetrics,fileNames,sleepStatus,phaseOfStudy,subjectStatus);
clear eegComputationalMetrics fileNames

load('Sleep2_Pre_results.mat')
[allUCBCasesResults ] = get_cvs_from_computational_metrics_mat(allUCBCasesResults,eegComputationalMetrics,fileNames,sleepStatus,phaseOfStudy,subjectStatus);
clear eegComputationalMetrics fileNames

phaseOfStudy = 'Post';
load('Sleep1_Post_results.mat')
[allUCBCasesResults ] = get_cvs_from_computational_metrics_mat(allUCBCasesResults,eegComputationalMetrics,fileNames,sleepStatus,phaseOfStudy,subjectStatus);
clear eegComputationalMetrics fileNames

load('Sleep2_Post_results.mat')
[allUCBCasesResults ] = get_cvs_from_computational_metrics_mat(allUCBCasesResults,eegComputationalMetrics,fileNames,sleepStatus,phaseOfStudy,subjectStatus);
clear eegComputationalMetrics fileNames


sleepStatus ='Awake';
phaseOfStudy = 'Pre';
load('Wake1_Pre_results.mat')
[allUCBCasesResults ] = get_cvs_from_computational_metrics_mat(allUCBCasesResults,eegComputationalMetrics,fileNames,sleepStatus,phaseOfStudy,subjectStatus);
clear eegComputationalMetrics fileNames

load('Wake2_Pre_results.mat')
[allUCBCasesResults ] = get_cvs_from_computational_metrics_mat(allUCBCasesResults,eegComputationalMetrics,fileNames,sleepStatus,phaseOfStudy,subjectStatus);
clear eegComputationalMetrics fileNames

phaseOfStudy = 'Post';
load('Wake1_Post_results.mat')
[allUCBCasesResults ] = get_cvs_from_computational_metrics_mat(allUCBCasesResults,eegComputationalMetrics,fileNames,sleepStatus,phaseOfStudy,subjectStatus);
clear eegComputationalMetrics fileNames

load('Wake2_Post_results.mat')
[allUCBCasesResults ] = get_cvs_from_computational_metrics_mat(allUCBCasesResults,eegComputationalMetrics,fileNames,sleepStatus,phaseOfStudy,subjectStatus);
clear eegComputationalMetrics fileNames




% Save the table to a CSV file
writetable(allUCBCasesResults, 'allPatinetsMetrics_Cz.csv');





