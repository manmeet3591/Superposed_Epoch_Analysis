clear all
clc
load('sea.csv')
load('ci.95.lower.csv')
load('ci.99.lower.csv')
load('ci.95.upper.csv')
load('ci.99.upper.csv')
load('lag.csv')
bar(lag,sea)
hold on
plot(lag,ci_95_lower)
hold on
plot(lag,ci_95_upper)
hold on
plot(lag,ci_99_lower)
hold on
plot(lag,ci_99_upper)
grid on