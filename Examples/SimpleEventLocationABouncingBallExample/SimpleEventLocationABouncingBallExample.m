%% Simple Event Location: A Bouncing Ball
% This example shows how to write a simple event function for use with an
% ODE solver. The example file |ballode| models the motion of a bouncing
% ball. The events function halts the integration each time the ball
% bounces, and the integration then restarts with new initial conditions.
% As the ball bounces, the integration stops and restarts several times.

%%
% The equations for the bouncing ball are
%
% $$\begin{array}{cl} y'_1 &= y_2\\ y'_2 &= -9.8 . \end{array}$$
%

%%
% A ball bounce occurs when the height of the ball $y_1(t)$ is equal to
% zero after decreasing. An events function that codes this behavior is
%
% <include>bounceEvents.m</include>
%

%%
% Type |ballode| to run the example and illustrate the use of an events
% function to simulate the bouncing of a ball.
ballode