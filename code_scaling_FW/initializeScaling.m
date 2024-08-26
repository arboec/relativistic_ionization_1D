function [v,width,u] = initializeScaling(SCALING)

if (SCALING == 1) %if scaling is ON
    %v = 0.625; % speed of mesh expansion. Used to calculate R(v,t)
    v = 0.58;%0.5 was defaul in July
    width = 0.5324; % with this width exp (which is exp(-(widt*x)^2)) coincides with my function at 1/2 level
    u = 2 * width*width; % parameter
else
    v = 0.0;
    width = 0.5324;
    u = 2 * width*width; % parameter
end

end

