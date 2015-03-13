This repository is no longer maintained.

Pyticles is an n-body physics code written in python, with a focus on the application of particle dynamics for continuum mechanics - SPH, DPD and related methods. The design prioritises flexibility over speed. Having said that, Cython provides a simple means to massive performance gains. Once features are implemented and tested I will be writing Cython versions of their most expensive inner loops.

The purpose of the code is to prototype algorithms and allow for interactively exploring models. More well known interactions such as springs and gravitation are also implemented.

Current execution speed is very slow - this is not so much a function of the inter-particle calculations as the sph rendering.

Some aspects of this development are related to my PhD work on smooth particle methods, so do get in touch if you want to contribute or just register your interest.

Acknowledgements:
Daniel Price' [SPLASH](http://www.astro.ex.ac.uk/people/dprice/pubs/index.html) code was a source of algorithms and inspiration.

Andrew Straw's [http://code.astraw.com/projects/motmot/wiki/pygarrayimage pygarrayimage module is used to connect numpy to Pyglet

Alex Holkner's [Pyglet](http://www.pyglet.org/) is used for the screen interface framework.

[Scipy](http://www.scipy.org/) and Numpy have provided the numerical nuts and bolts that make doing this in python (without writing my own array code... ) even remotely possible.

Peter Daivis' molecular dynamics code was a source for various algorithms.

The Lepton high performance particle engine has been a source of ideas and inspiration http://code.google.com/p/py-lepton/.

The pronunciation is 'Pie-tic-lees', as in Patrocles, Heracles