# FortSEIR
FortSEIR is a Fortran based library to simulate disease spread amongst different demographics in a population. It makes use of Intel's MKL library for fast matrix operations and uses a 4-th order, error controlled, Runge-Kutta method to solve the coupled non-linear ODEs that describe the spread of disease in the population. More details about the equations can be found in the attached file SEIR_age_classified.pdf

The code is free to use and is still under development. The final version will be notified soon. There are two versions: Age-structured and non-age-structured (I call this toy).

# Data for COVID19 India (https://www.kaggle.com/sudalairajkumar/covid19-in-india/data)
Data for COVID19 outbreak in India is in covid_Data_India (upto April 8th, 2020). Data updated daily at link.
Import data from covid_19_india.xlsx as column vectors in Matlab through uiimport and save the file with original column names. Use read_COVID_Data.m in Matlab_codes to read data.
