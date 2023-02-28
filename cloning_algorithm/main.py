#!/usr/bin/env python

import numpy as np
import multiprocess as multiprocessing
import copy
import h5py


class cloning_algorithm():
    '''
    This class implements a cloning algorithm similar to those described in 

    [1] "Experimental Measurement of Relative Path Probabilities and 
        Stochastic Actions", J. Gladrow, U. F. Keyser, R. Adhikari, and 
        J. Kappler.
        Physical Review X, vol. 11, p. 031022, 2021. 
        DOI: 10.1103/PhysRevX.11.031022
    [2] "Measurement of irreversibility and entropy production via the 
        tubular ensemble", J. Kappler and R. Adhikari. 
        Physical Review E, vol. 105 (4), p. 044107, 2022. 
        DOI: 10.1103/PhysRevE.105.044107

    ##########################################################
    # Problem description and short summary of the algorithm #
    ##########################################################

    The cloning algorithm allows to measure probabilities of rare events of 
    Markov processes. The cloning algorithm can be applied to rare events that
    can be formulated as survival probabilities.

    For example, let us consider the following question: For the Wiener 
    process starting at x0 = 0, what is the survival probability P to never 
    leave the interval [-1,1] until time T?

    One would naively measure this survival probability by generating a large
    number of time series (in our example: realizations of the Wiener 
    process), and evaluating the fraction of time series that fulfill the
    survival condition up to the final time (in our example: count how many 
    time series have never left the interval [-1,1] and divide by the total
    number of time series that were generated).
    
    The problem with this naive approach is that the survival probability 
    often decays exponentially. In our example it decays as as -log(P) ~ T,
    so that one needs an unfeasibly large number of time series to get an
    acceptable estimate of the survival probability.

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

    Note that the paragraphs above are a slightly simplified description. In
    particular, the actual cloning algorithm estimates a reasonable value for
    n_clones at each iteration based on the recent decay of the survival 
    probabiliy.

    For a more detailed explanation for the cloning algorithm can be found
    in Fig. 2 and Appendix B of Ref. [1], as well as Fig. 1 and Appendix B 
    of Ref. [2]. 
    
    ##############################################
    # How does the algorithm get the time series #
    # and evaluate the survival condition?       #
    ##############################################

    An instance of the cloning_algorithm takes as input instances of two 
    user-defined classes:
      1. A time_series_generator that generates stochastic sample time series
      2. A condition_detector, that checks up to which timestep a sample time 
         series fulfills a condition
    These have to take input and output in the following form (if an existing
    simulation code has a different way of taking inputs/outputting results,
    one can write a simple wrapper to put the results into the format 
    described below):

    1. A time_series_generator, which is an instance of any class that 
       contains a method 
            x = time_series_generator.get_time_series(x0,t0)
       that takes as input an array x0 of length n_dim (with n_dim the 
       dimensionality of the stochastic trajectory), and a float t0. The 
       tuple (x0,t0) represents the initial condition and initial time of a 
       stochastic trajectory. 
       The method returns an array x of shape (n_steps+1,n_dim), with n_steps
       the number of steps per iteration of the cloning algorithm. 
    2. A condition_detector, which is an instance of any class that contains
       the following two methods:
       - leaving_index = condition_detector.evaluate_condition(x,t):
         This method takes as input an array x of shape (n_steps+1,n_dim) with
         a trajectory and corresponding times t of shape (n_steps+1).
         It then returns an integer leaving_index such that 
         (x[leaving_index], t[leaving_index]) is the first datapoint at which
         the desired condition is not fulfilled. If leaving_index = -1, then 
         the condition is fulfilled until the final step of the time series.
       - observables = condition_detector.evaluate_observables(results):
         This method takes as input a dictionary which contains:
         * results['t'] = array of shape (n_steps + 1) with times,
         * results['x'] = array of shape (n_clones, n_steps + 1, n_dim) with 
           time series,
         * results['leaving_index'] = array of shape (n_clones), where the 
           i-th element is the index at which trajectory results['x'][i] does
           not fulfill the desired condition for the first time,
         * results['final_position'] = array of shape (n_clones, n_dim) where
           the i-th element is final position of trajectory i
         The function can use any of those inputs to calculate any observables
         one might be interested in, and returns those observables in a 
         dictionary observables.

    '''

    def __init__(self,time_series_generator,condition_detector,parameters):
        #
        self.set_standard_parameters()
        #
        self.set_parameters(parameters=parameters)
        #
        self.set_time_series_generator(time_series_generator=\
                                            time_series_generator)
        #
        self.set_condition_detector(condition_detector=condition_detector)
        #
        self.create_return_data_dictionary()
        #

    def set_standard_parameters(self):
        #
        self.n_clones = 10000
        self.n_clones_max = 1000000
        self.n_chunksize = 1000
        self.n_cores = 1
        self.n_target = 1000
        self.n_fit = 50
        self.verbose = False
        self.dt = 1.

    def set_condition_detector(self,condition_detector):
        self.condition_detector = condition_detector

    def set_time_series_generator(self,time_series_generator):
        self.time_series_generator = time_series_generator

    def set_parameters(self,parameters):
        #
        try:
            self.n_clones = parameters['n_clones']
        except KeyError:
            pass
        #
        try:
            self.n_clones_max = parameters['n_clones_max']
        except KeyError:
            pass
        #
        try:
            self.n_target = parameters['n_target']
        except KeyError:
            pass
        #
        try:
            self.n_fit = parameters['n_fit']
        except KeyError:
            pass
        #
        try:
            self.n_iterations = parameters['n_iterations']
        except KeyError:
            pass
        #
        try:
            self.n_chunksize = parameters['n_chunksize']
        except KeyError:
            pass
        #
        try:
            self.n_cores = parameters['n_cores']
        except KeyError:
            pass
        #
        try:
            self.n_steps = parameters['n_steps']
        except KeyError:
            pass
        #
        try:
            self.initial_conditions = parameters['initial_conditions']
        except KeyError:
            pass
        #
        try:
            self.t0 = parameters['t0']
            self.t = self.t0
        except KeyError:
            pass
        #
        try:
            self.dt = parameters['dt']
        except KeyError:
            pass
        #
        try:
            self.verbose = parameters['verbose']
        except KeyError:
            pass
        #
        self.t = np.arange(self.n_steps+1,dtype=float)*self.dt + self.t0


    def get_parameters(self):
        parameters = {'n_clones':self.n_clones,
                    'n_clones_max':self.n_clones_max,
                    'n_iterations':self.n_iterations,
                    'n_chunksize':self.n_chunksize,
                    'n_cores':self.n_cores,
                    'n_steps':self.n_steps,
                    'n_target':self.n_target,
                    'n_fit':self.n_fit,
                    't0':self.t0,
                    'dt':self.dt
                    }
        return parameters

    def create_return_data_dictionary(self):
        #
        self.return_data = {'t':np.arange(self.n_iterations*self.n_steps+1)\
                                                                *self.dt,
                            'remainers':np.array([self.n_clones],dtype=int),
                            'survival_probability':np.array([1.],dtype=float),
                            'exit_rate':np.array([],dtype=float),
                            'n_clones':[],
                            'final_positions':[],
                            'observables':{}
                            }
        self.return_data['parameters'] = self.get_parameters()

    def get_return_data_dictionary(self):
        return self.return_data

    def get_time_series_and_evaluate_condition(self,x0):
        #
        x = self.time_series_generator.get_time_series(x0=x0,t0=self.t[0])
        #
        leaving_index = self.condition_detector.evaluate_condition(x=x,
                                                                t=self.t)
        #
        results = {'x':x.copy(),
                   'leaving_index':leaving_index,
                   'final_position':x[leaving_index]}
        #
        return results

    def append_remainers(self,results):
        #
        leaving_indices = results['leaving_index']
        #
        bin_edges = np.arange(0,self.n_steps+1)+0.5
        # get number of remaining trajectories as function of timestep
        hist, bin_edges = np.histogram(np.array(leaving_indices),
                                        bins=bin_edges)
        # hist_leaving_indices[0] == 0, since all trajectories initially
        # fulfilled the condition. We thus start with index 1 for the time
        # series of the remainers:
        remainers = self.n_clones - np.cumsum(hist)
        remainers_normalized = remainers.astype(float)/self.n_clones
        #
        # len(remainers) = self.n_steps
        self.return_data['remainers'] = np.concatenate(
                        (
                            self.return_data['remainers'],
                            remainers
                        ))
        self.return_data['survival_probability'] = np.concatenate(
                        (
            self.return_data['survival_probability'],
            remainers_normalized*self.return_data['survival_probability'][-1]
                        ))
        #
        return 0

    def append_exit_rate(self):
        #
        remainers = self.return_data['survival_probability']
        #
        if len(remainers) == self.n_steps + 1: # means we are at the end
            # of the first iteration
            current_rate = np.zeros(self.n_steps,dtype=float)
            #
            current_rate[1:] = remainers[2:]-remainers[:-2]
            current_rate[1:] /= -2*remainers[1:-1]
            #
            current_rate[0] = remainers[1]-remainers[0]
            current_rate[0] /= -remainers[0]
            #
        else:
            i = self.n_steps+2
            recent_remainers = remainers[-i:]
            current_rate = recent_remainers[2:]-recent_remainers[:-2]
            current_rate /= -2*recent_remainers[1:-1]
        #
        self.return_data['exit_rate'] = np.concatenate(
                        (
                            self.return_data['exit_rate'],
                            current_rate
                        ))
        
    def calculate_final_exit_rate(self):
        #
        remainers = self.return_data['survival_probability']
        #
        final_value = -(remainers[-1] - remainers[-2])/remainers[-1]
        #
        self.return_data['exit_rate'] = np.concatenate(
                        (
                            self.return_data['exit_rate'],
                            np.array([final_value])
                        ))
        #
        self.return_data['exit_rate'] = self.return_data['exit_rate']/self.dt
        #

    def append_observables(self,results):
        #
        observables = self.condition_detector.evaluate_observables(results)
        #
        for key, value in observables.items():
            try: 
                self.return_data['observables'][key] = np.concatenate((
                        self.return_data['observables'][key],
                        value
                ),axis=0)
            except KeyError:
                self.return_data['observables'][key] = value

    def update_n_clones(self):
        #
        self.return_data['n_clones'].append(self.n_clones)
        #
        # estimate decay in next iteration from recent exit rate 
        exit_rate = self.return_data['exit_rate']
        if len(exit_rate) < self.n_fit:
            recent_exit_rate = exit_rate
        else:
            recent_exit_rate = exit_rate[-self.n_fit:]
        #
        n_recent_exit_rate = len(recent_exit_rate)
        x_fit = np.arange(n_recent_exit_rate,dtype=float) \
                            - (n_recent_exit_rate-1)
        lin_fit = np.polyfit(x_fit,recent_exit_rate,deg=1)
        #
        if lin_fit[0] < 0:
            self.n_clones = int( self.n_target \
                                * np.exp(lin_fit[1]*self.n_steps) 
                                )
        else:
            self.n_clones = int( self.n_target  \
                                 * np.exp(lin_fit[1]*self.n_steps \
                                 + lin_fit[0]/2*(self.n_steps)**2) 
                                 )
        #
        return 0
    
    def check_size_of_n_clones(self):
        #
        if self.n_clones > self.n_clones_max:
            #
            warning_msg = ("Warning: n_clones = {0} exceeds "
                    "n_clones_max = {1}. Setting n_clones = "
                    "n_clones_max = {1}")
            print(warning_msg.format(self.n_clones,
                                        self.n_clones_max))
            #
            self.n_clones = self.n_clones_max
        return 0

    def update_initial_conditions(self,results):
        #
        # append old initial conditions to return dictionary
        self.return_data['final_positions'].append(
                                    copy.deepcopy(self.initial_conditions)
                                                )
        #
        # collect final positions of those trajectories 
        # that remained in the tube
        mask = (results['leaving_index'] < 0)
        self.initial_conditions = results['final_position'][mask]

    def get_number_of_survivors(self,results):
        #
        leaving_index = results['leaving_index']
        return np.sum((leaving_index < 0 ))

    def propagate(self):
        #
        # draw initial conditions for the simulations
        indices = np.random.randint(0,len(self.initial_conditions),
                                        self.n_clones)
        #
        # Run simulations. Multithreading is implemented based on this post:
        # https://stackoverflow.com/a/9786225
        pool = multiprocessing.Pool(self.n_cores)
        #
        list_of_results = list(
                    pool.imap(self.get_time_series_and_evaluate_condition,
                            [self.initial_conditions[e] for e in indices],
                            chunksize=self.n_chunksize)
                            )
        pool.close()
        pool.join()
        #
        # turn list of dictionary into dictionary of arrays
        # results = dictionary that contains both the time series of the
        #           current iteration, as well as the information up to which
        #           timestep each time series fulfills the condition
        for i,e in enumerate(list_of_results):
            if i == 0:
                results = {key:[value] for key, value in e.items() }
            else:
                for key, value in e.items():
                    results[key].append(value)
        results['t'] = (self.t).copy()
        #
        for i,e in results.items():
            results[i] = np.array(e)
        #
        results['n_survivors'] = self.get_number_of_survivors(results)
        #
        if results['n_survivors'] < 1:
            return results
        #
        self.append_remainers(results)
        #
        self.append_exit_rate()
        #
        self.append_observables(results)
        #
        self.update_initial_conditions(results)
        #
        self.update_n_clones()
        #
        self.check_size_of_n_clones()
        #
        return results

    def run(self,filename=None):
        #
        status_msg = ("Running iteration {0} of {1} with {2} clones")
        finish_msg = ("Finished all {0} interations. Mean number of clones "
                        "used: {1}.                 ")
        #
        for i in range(self.n_iterations):
            if self.verbose:
                print(status_msg.format(i+1,self.n_iterations,self.n_clones),
                            end='\r')
            #
            #
            results = {'n_survivors':0}
            #
            #
            while (results['n_survivors'] == 0):
                results = self.propagate()
                if results['n_survivors'] == 0:
                    # 
                    self.n_clones *= 2
                    #
                    warning_msg = ("Warning: In iteration {0} of {1}, no"
                            " time series fulfilled the condition until the"
                            " final time of the iteration. Rerunning the"
                            " iteration with n_clones = {2}.")
                    print(warning_msg.format(i+1,
                                            self.n_iterations,
                                            self.n_clones))
                    #
                    self.check_size_of_n_clones()
            #
            self.t += self.n_steps*self.dt
        #
        self.calculate_final_exit_rate()
        #
        if self.verbose:
            print(finish_msg.format(self.n_iterations,
                                np.mean(self.return_data['n_clones'])),
                    end='\n')
        #
        if filename != None:
            self.save(filename=filename)
        #
        return self.return_data

    def save(self,filename='cloning_results.h5'):
        #
        with h5py.File(filename, 'w') as hf:
            for key, value in self.return_data.items():
                if key.lower() in ['final_positions']:
                    # see https://stackoverflow.com/a/37218874
                    grp = hf.create_group(key)
                    #value = np.array(value,dtype=object)
                    for i,e in enumerate(value):
                        grp.create_dataset(str(i),data=np.array(e))
                elif key.lower() in ['parameters','observables']:
                    #
                    grp = hf.create_group(key)
                    #
                    for i,e in value.items():
                        grp.create_dataset(i,data=e)
                else:
                    hf.create_dataset(key, 
                                    data=value, 
                                        )



