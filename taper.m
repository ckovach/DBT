classdef taper

% The taper class defines some properties of window tapers. 
%
%Properties:
%
% generator:  Handle defining a function, h, used to generate the taper.
%             h is defined over [0 1] and typically has h(0) = 1 and h(1) = 0. 
%             Default is @(x)1/2*(1+cos(x*pi)). This is not the same as
%             the taper function,
%             (see the following properties)
%
% summation_order:  Require that the taper, f(x), satisfy 
%                   f(x)^p + f(1-x)^p = 1.  Default is p = 2. This creates
%                   a function  g(x) = h(x) for x>=0 and g(x) = (1-h(1-abs(x))^p)^1/p for 
%                   x < 0. For p = 0 g(x) = h(x).
%
% symmetric:   Require that the taper be symmetric if true. If true,
%              then the taper function is
%                   f(x) = (g(x)+g(-x))/2.
%              Default is true.
%Methods:
%
% obj = taper; %Returns taper with default properties
%
% val = obj.make(x).  Returns taper values at window points in the
%                     range [-1 1]. 
% 
%
    
% C Kovach 2013

    properties
        
        
        generator = @(x)cos(x*pi/2);
        
        summation_order=2;
        
        symmetric = false;
    end
    
    methods
        
     
        
        function [t,f] = make(me,x)
           
            h = me.generator;
            p = me.summation_order;
            if p > 0
              %  g = @(x) (x>=0).*h(abs(x)) + (x<0).*(1-(h(1-abs(x)).^p)).^(1/p);
                g = @(x) (mod(x,2)-1<0 | x==0).*h(mod(x,2)) + (mod(x,2)-1>=0 & x~=0).*(1-h(mod(x,2)-1).^p).^(1/p);
                               
            else
                g = h;
            end
            
            if me.symmetric
                f = @(x) ((g(x).^p + g(-x).^p)/2).^(1/p);
            else
                f = g;
            end
               
            t = f(x);
            
        end
    end
end
