function [txt, y, res, sigcor, fitobject] = mkfit(t ,sig, BleachModel)
%MKFIT is fitting data according to different models to estimate the
%bleaching occuring during fluorescent time lapse experiments.
% 
%   MKFIT(t,sig,BleachModel) is fitting sig(t) using different models
%   specified by BleachModel. t and sig have to be vectors of equal
%   dimension. BleachModel can be: none, linear, single exponential or
%   double exponential.
%
%   [txt,y,res,sigcor,fitobject] = MKFIT(...) where tyt is a string
%   providing the parameters of the fitting model, y is a vector of the
%   fitted data from sig, res is the residue, sigcor is the signal
%   corrected for the bleach. Finally  fitobject is fit object object given
%   by the model and the data.
%
%
%   Nicolas Liaudet
%   Bioimaging Core Facility - UNIGE
%   https://www.unige.ch/medecine/bioimaging/en/bioimaging-core-facility/
% 
%   v1.0 05-Jul-2018 NL


t = t';
sig = sig';
switch BleachModel
    case 'none'
        fitobject = 'none';
        y = sig;
        sigcor = sig;
        txt = ['s = s'];
    case 'linear'
        fitobject = fit(t,sig,'poly1');
        y = fitobject.p1*t+fitobject.p2;
        sigcor = sig-y+y(1);
        txt = ['s = ',...
            num2str(fitobject.p1) 't+' num2str(fitobject.p2)];
        
    case 'single exponential'
        fitobject = fit(t,sig,'exp1');
        y = fitobject.a*exp(fitobject.b*t);
        sigcor = sig-y+y(1);
        txt = ['s = ',...
            num2str(fitobject.a) 'e^{' num2str(fitobject.b) 't}'];
        
    case 'double exponential'
        fitobject = fit(t,sig,'exp2');
        y = fitobject.a*exp(fitobject.b*t)+...
            fitobject.c*exp(fitobject.d*t);
        sigcor = sig-y+y(1);
        txt = ['s = ',...
            num2str(fitobject.a) 'e^{' num2str(fitobject.b) 't} + ',...
            num2str(fitobject.c) 'e^{' num2str(fitobject.d) 't}'];
end
res = sig-y;
end
