# Cloning algorithm for measuring rare events from stochastic time series

This module implements a cloning algorithm similar to those described in Refs. <a href="#ref_1">[1]</a>, 
<a href="#ref_2">[2]</a>.

The cloning algorithm allows to <b>measure probabilities of a class of rare events of 
Markov processes</b>. More explicitly, the cloning algorithm can be used to measure probabilities of rare events that
can be formulated as survival probabilities.

## <a id="problem_statement"> Problem: It is often infeasible to straightforwardly measure probabilities of rare events

For example, consider the following question: For a freely diffusion particle (Wiener 
process) starting at $X_{t=0} = 0$, what is the survival probability $P(T)$ to never 
leave the interval $[-1,1]$ until some large time $t = T$?

One would naively measure this survival probability by generating a large
number of sample time series (realizations of the stochastic process), 
and evaluating the fraction of time series that fulfill the
survival condition. 
In our example this means first generating $N$ realizations of the Wiener process (each of which starts at the origin),
and then counting how many 
sample time series of duration $T$ have never left the interval
$[-1,1]$. If of the $N$ generated realizations,  $M$ have never the interval, the ratio $M/N$ is a reasonable estimate
for the survival probability, $P(T) \approx M/N$.

The problem with this approach is that the survival probability 
often decays exponentially. In our example it decays as $\log(P) \sim -T$,
so that one needs an exponentially large number of sample time series to get an
acceptable estimate of the survival probability. It is often infeasible to generate
and analyze such a large number of time series.

## <a id="cloning_algorithm"> Solution: The cloning algorithm

The cloning algorithm is a means for overcoming the exponential decay of the number of 
surviving samples.

To do this, the algorithm divides the time interval $I = [0,T]$ into
 $n_{\mathrm{iterations}}$ shorter time intervals of duration $\Delta T := T/n_{\mathrm{iterations}}$,
 and treats the survival probability for each interval separately:
 
1. For the first interval $I_0 = [0, \Delta T]$, $N^{(0)}_{\mathrm{clones}}$ time series are generated using the
 initial conditions of the problem. Those time series are then used to estimate the survival probability naively (as discussed in the previous section) for the duration $I_0$.
2. For the second time interval $I_1 = [\Delta T, 2 \Delta T]$,  $N^{(1)}_{\mathrm{clones}}$ time series are then generated,
with initial conditions drawn from the final positions of the survivors from the previous iteration.
The resulting time series are then used to estimate the survival probability in the time interval $I_1$.
4. Step 2 is repeated until the survival probability has been obtained up to the final time $T$.

The "cloning" in the name of the  algorithm comes from the fact that at the end of
each iteration, we take the surviving final states and "clone" them by drawing
initial conditions for the next iteration. Note that to be able to use the final positions of
each batch of trajectories as initial conditions for the next batch, we assume
that the time series is Markovian.

By considering the survival problem on each interval $I_i$ separately, and periodically increasing
the number of sample trajectories, we overcome the 
exponential temporal decay of sample time series that fulfill the 
survival condition. 

For chosen duration $\Delta T$, the implemented cloning algorithm 
uses the recent decay of the survival probabiliy to estimates a reasonable value for
the number of samples
$N^{(i)}_{\mathrm{clones}}$ at each iteration.

For a more detailed explanation for the cloning algorithm can be found
in Fig. 2 and Appendix B of <a href="#ref_1">Ref. [1]</a>, as well as Fig. 1 and Appendix B 
of Ref. <a href="#ref_2">Ref. [2]</a>. 


## <a id="examples"> Examples
  
### <a id="example_sojourn"> Sojourn probability

Example notebooks are provided in the folder [examples/](examples):

* **sojourn probability.ipynb**: In this notebook, we measure the survivial probability for freely diffusive stochastic dynamics to remain within a tube of radius $R$ around a given reference path $\varphi(t)$ (in this context the surviving probability is also called sojourn probability). We evaluate the decay rate corresponding to the survival probability, and compare the measured decay rate to theoretical predictions based on the Onsager-Machlup stochastic action.

 More notebooks to be added.


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
