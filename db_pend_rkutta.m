clc;
clear all;



stv_in=[4*pi/180,0,0,0]
ts=linspace(0,2,100);

[tout,yout]=ode45(@dydt_dbp,ts,stv_in);
