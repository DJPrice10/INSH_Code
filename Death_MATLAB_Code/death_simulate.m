
function y_ret = simulate_death(b1,t,N)
%%
% PLEASE NOTE THAT IT CAN BE VERY HELPFUL TO WRITE THE CODE TO SIMULATE
% FROM THE MODEL IN A LOWER LEVEL LANGUAGE LIKE C AND CALL IT INTO MATLAB.
% TO ILLUSTRATE THE METHODOLOGY HERE WE HAVE INCLUDED MATLAB CODE ONLY.
%
% simulate data from the death model
%
% input:
%   t - times to record epidemic
%   b1 - value of parameter
%   N - initial number of susceptibles
%
% output:
%   y_ret - number of infecteds at times t
%%
y_ret = zeros(length(t),1);
t_curr = 0;
obs_counter=1;
S_t = N;
I_t = 0;
num_obs = length(t);

while(1)
    % simulate time until next event
    t_curr = t_curr + exprnd(1/(b1*S_t));
    while (obs_counter <= num_obs && t_curr > t(obs_counter))
        y_ret(obs_counter) = I_t;
        obs_counter = obs_counter+1;
    end
    if (t_curr > t(num_obs))
        return;
    end
    % now update the states
    I_t = I_t+1; S_t = S_t-1;
    if (S_t==0)
        % then epidemic is OVER, fill in the rest of observations
        for i = obs_counter:num_obs
            y_ret(i) = N;
        end
        return;
    end
end
	

end