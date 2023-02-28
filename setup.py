from distutils.core import setup


description = ('python implementation of a cloning algorithm'
                'for measuring the probability of rare events')
setup(
    name='cloning_algorithm',
    version='0.0.1',
    url='https://github.com/juliankappler/',
    author='Julian Kappler',
    license='MIT',
    description=description,
    long_description=description,
    platforms='tested on macOS, should work on all platforms',
    libraries=[],
    packages=['cloning_algorithm'],
)
