% ======================================================================
%> @file main.m
%> @brief Main script for heterogeneous mechanical test design
%> @details
%>
%> @author Mafalda Gonçalves
%> @date since 2020
% ======================================================================
clear; close all; clc
tic

inputfile_name = 'Input/Case1_total_homogeneous.dat';
[input_data] = inputfile_processing(inputfile_name);

% %Main function
boosted_TopOpt(input_data); 
toc


%======================================================================
%>@file main.m
%>@brief Main script for heterogeneous mechanical test design
%>@details
%>  
%>@author Mafalda Gonçalves
%>@date since 2020
%======================================================================