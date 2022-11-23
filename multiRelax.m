function [T,X] = multiRelax(x0, totLenStimuli)
    % This function simulate the networok without stimuli intervention

	X = zeros((totLenStimuli*10)+1,length(x0));
    T = zeros((totLenStimuli*10)+1,1);
    t0 = 0;
    h = 0.01;
    x1 = x0;
    k = 1;
    for i = 1:totLenStimuli
        tf = t0+h*10;
        tspan = t0:h:tf;
        [t,x]=ode23tb(@f,tspan,x1);
        x1 = x(end,:).';
        X(k:k+10,:) = x;
        T(k:k+10,:) = t;
        k = k+10;
        t0 = tf;
    end
end