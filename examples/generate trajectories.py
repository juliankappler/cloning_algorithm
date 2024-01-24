#!/usr/bin/env python

import numpy as np
import h5py
from pathlib import Path # for creatin the directory where the sample time 
# series are saved

class stoch_ode:
    '''
    This class generates stochastic time series for given number of timesteps 
    and initial state
    '''

    def __init__(self,parameters):
        #
        self.n_steps = parameters['n_steps']
        self.dt = parameters['dt']
        self.D0 = parameters['D0']
        self.a0 = parameters['a0']
        self.a1 = parameters['a1']
        #
        try:
            self.L = parameters['L']
        except KeyError:
            self.L = 1.
        #
        self.t = np.arange(self.n_steps+1)*self.dt
        #
        self.a = lambda x: self.a1/self.L \
                            + 4*self.a0*x/self.L**2 \
                            - 4*self.a0*x**3/self.L**4
        self.D = lambda x: self.D0*( 5 - 1*np.cos(np.pi*x/self.L))/4.
        #

    def simulate(self,x0,
                    t0=0,
                    filename=None):
        #
        traj = np.zeros([self.n_steps+1,1],
                        dtype=float)
        traj[0] = x0
        #
        self.rng = np.random.default_rng()
        noise = self.rng.normal(size=self.n_steps)
        #
        for i,cur_x in enumerate(traj[:-1]):
            traj[i+1] = cur_x \
                        + self.a(cur_x) * dt \
                        + np.sqrt(2*self.D(cur_x)*self.dt) * noise[i] 
        #
        t = t0 + self.t
        #
        result = {'t':t,'x':traj}
        #
        if filename is not None:
            self.save(dictionary=result,
                        filename=filename)
        #
        return result

    def save(self,
                dictionary,
                filename):
        #
        with h5py.File(filename, 'w') as hf:
            for key, value in dictionary.items():
                hf.create_dataset(key, 
                                data=value, 
                                    )




# parameters for simulation
t0 = 0.
dt = 5e-5
D0 = 1.
a0 = 2.
a1 = a0/20.
n_steps = 500000

N_simulations = 20

output_directory = "./sample_data/"

parameters = {'D0':D0,
                'dt':dt,
                'n_steps':n_steps,
                'a0':a0,
                'a1':a1,
                }

simulator = stoch_ode(parameters)

if __name__ == '__main__':
    #
    Path(output_directory).mkdir(parents=True, exist_ok=True)
    #
    rng = np.random.default_rng()
    x0s = 2*(rng.random(size=N_simulations)-0.5)
    #
    for i in range(N_simulations):
        print('Running simulation {0} of {1}     '.format(i+1,N_simulations),
                        end='\r')
        #
        filename = '{0}/trajectory_{1:05d}.h5'.format(output_directory,i)
        results = simulator.simulate(x0=x0s[i],
                            filename=filename)