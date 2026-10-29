# RTMB-PSC-2026

This repository contains materials for the RTMB workshop hosted by the Pacific Salmon Commission in November 2026.

This is scheduled as a two-day workshop from 9am-4:30pm. The schedule will provide plenty of time for breaks and lunch. 
The location will be the Pacific Salmon Commission Office 600-1155 Robson St, Vancouver BC V6E 1B5.

To prepare for the workshop:

 - Complete module_0 pre-workshop assignment
 - Install RTMB
 - Download the materials provided in this repository (the day before the course as it will be updated before then)

All materials for the workshop will be in this GitHub repository. If you're familiar with Git/GitHub, you already know how to get all the materials on your computer. 
If you're not, go [here](https://github.com/Pacific-salmon-assess/RTMBWorkshop_Nov2026/), click the (green) "Code" button, and choose the "Download ZIP" option.

## Background for the workshop

This workshop will focus on the writing stock assessment models using the R package `RTMB`. The material assumes attendees have basic knowledge of 
stock assessment models, comfortable with the program R, and can create a working model in code based on a likelihood equation.

## Tentative Schedule

Day 1: Monday November 16

1. (9:00 – 9:50) Welcome and course introduction, introduction to maximum likelihood estimation, and assignment review.
2. (10:00 - 10:50) Writing and fitting a stock recruit likelihood in R. Improve that likelihood with RTMB and compare. Compare with standard R packages. 
3. (11:00 – 11:50) Introduction to RTMB and automatic differentiation (AD). Compare AD with finite differencing. What is a gradient and a Hessian. 
4. (13:00 – 13:50) Add random effects to your likelihood. Laplace and RTMB 
5. (14:00 – 14:50) Practice writing likelihoods with random effects and fit them. Hierarchical stock recruitment models.  Compare output with `glmmTMB`.
6. (15:00 – 16:30) Additional random effect structures (e.g. time varying productivity). 

Day 2: Tuesday November 17

1. (9:00 – 9:50) Variance estimation for maximum likelihood. Delta method to approximate parameters. Using RTMB to report variance and confidence intervals. 
2. (10:00 – 10:50) RTMB functionality. Learn little details about how to use RTMB effectively. 
3. (11:00 – 11:50) Writing helper functions and atomics (e.g. Lambert W). Best coding practices. 
4. (13:00 – 16:30) Practice writing complex likelihood functions in RTMB 
