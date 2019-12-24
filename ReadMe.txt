This code estimates SV models through Bayesian estimation methods
with posterior distribution for the parameters approximated 
through MCMC simulations.


How to use the code:

First compile,

sh compile.sh

It creates 8 executables:
SVEstim
SVJEstim
SV2FEstim
SVJ2FEstim
SVDUREstim
SV2FDUREstim
SVJDUREstim
SV2FDUREstim
SVJ2FDUREstim

How to run:

./executable returns_file num_iterations num_warmup num_thin output_file (without durations)
./executable returns_file durations_file num_iterations num_warmup num_thin output_file (with durations)

Examples:

./SVEStim y5aapl.txt 10000 50000 10 sv5aapl.txt
./SV2FDUREstim y5aapl.txt d5aapl.txt 10000 50000 10 sv2fdur5aapl.txt


Data files
Returns:

y5aapl.txt
y15aapl.txt
y5ge.txt
y15ge.txt

Durations:

d5aapl.txt
d15aapl.txt
d5ge.txt
d15ge.txt



