# Cloning algorithm for measuring rare events from stochastic time series

This module implements a cloning algorithm similar to those described in Refs. <a href="#ref_1">[1]</a>, 
<a href="#ref_2">[2]</a>.

The cloning algorithm allows to <b>measure probabilities of a class of rare events of 
Markov processes</b>. More explicitly, the cloning algorithm can be used to measure probabilities of rare events that
can be formulated as survival probabilities.

## <a id="problem_statement"> Problem: It is often infeasible to directly sample rare events

For example, consider the following question: For a freely diffusion particle (Wiener 
process) starting at <i>x(t=0) = 0</i>, what is the survival probability <i>P</i> to never 
leave the interval <i>[-1,1]</i> until some large time <i>T</i>?

One would naively measure this survival probability by generating a large
number of sample time series, and evaluating the fraction of time series that fulfill the
survival condition. 
In our example this means generating <i>N</i> sample trajectories of the Wiener process (each of which starts at the origin),
counting how many 
sample time series of duration <i>T</i> have never left the interval
<i>[-1,1]</i> (say <i>M</i>) and using the ratio <i>M/N</i> as an estimator 
for the survival probability, <i>P approx M/N</i>.

The problem with this naive approach is that the survival probability 
often decays exponentially. In our example it decays as <i>-log(P) ~ T</i>,
so that one needs an exponentially large number of time series to get an
acceptable estimate of the survival probability. It is often infeasible to generate
and analyze such a large number of time series.

## <a id="cloning_algorithm"> Solution: The cloning algorithm

The cloning algorithm is a means for overcoming the exponential decay of 
surviving samples.

The cloning algorithm solves this problem by dividing the time T into
shorter n_iterations shorter time intervals of duration T/n_iterations.
At the first time interval, n_clones time series are generated using the
given initial conditions, and used to estimate the survival probability
until time T/n_iterations. 
For the second time interval, n_clones time series are then generated
using as initial conditions that survivors from the previous iteration.
The resulting time series are then used to estimate the survival 
probability in the time interval [T/n_iterations, 2 * T/n_iterations].
This process is repeated until the survival probability has been obtained
up to time T.

The "cloning" in the algorithm name comes from the fact that at the end of
each iteration, we take the surviving final states (which are for 
well-chosen durations T/n_iterations about one order of magnitutde 
smaller in number than n_clones) and "clone" them by drawing n_clones
initial conditions for the next iteration. By this we overcome the 
(exponential) temporal decay of sample time series that fulfill the 
survival condition. That this works is of course due to the assumed Markov
property of the time series.

The paragraphs above are a slightly simplified description. In
particular, the actual cloning algorithm estimates a reasonable value for
n_clones at each iteration based on the recent decay of the survival 
probabiliy.

For a more detailed explanation for the cloning algorithm can be found
in Fig. 2 and Appendix B of <a href="#ref_1">Ref. [1]</a>, as well as Fig. 1 and Appendix B 
of Ref. <a href="#ref_2">Ref. [2]</a>. 


## <a id="examples"> Examples
  
### <a id="example_sojourn"> Sojourn probability

More details coming soon. For now, see the jupyter notebook [examples/sojourn probability.ipynb](examples/sojourn%20probability.ipynb).



## <a id="installation">  Installation

To install the module cloning_algorithm, as well as its requirements ([NumPy](https://numpy.org/), [multiprocess](https://pypi.org/project/multiprocess/), and [h5py](https://pypi.org/project/h5py/), clone this repository and run the installation script:

```bash
>> git clone https://github.com/juliankappler/cloning_algorithm.git
>> cd cloning_algorithm
>> pip install -r requirements.txt
>> pip install .
```


## <a id="references"> References

<a id="ref_1">[1] **Experimental Measurement of Relative Path Probabilities and Stochastic Actions**. J. Gladrow, U. F. Keyser, R. Adhikari, and J. Kappler. Physical Review X, vol. 11, p. 031022, 2021. DOI: [10.1103/PhysRevX.11.031022](https://doi.org/10.1103/PhysRevX.11.031022).</a>

<a id="ref_2">[2] **Measurement of irreversibility and entropy production via the tubular ensemble**. J. Kappler and R. Adhikari. Physical Review E, vol. 105 (4), p. 044107, 2022. DOI: [10.1103/PhysRevE.105.044107](https://doi.org/10.1103/PhysRevE.105.044107).</a>
