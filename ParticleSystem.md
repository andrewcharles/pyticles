ParticleSystem objects store the states of the particles. Currently their update() method is responsible for integrating the model forward in time.

One particle system can be operated on by multiple forces, and can have multiple interaction lists. I am currently investigating using a sorted verlet list to replace
some of the multiple neighbour list functionality.

An aim is to support coupling between different particle systems. For example you might have a simple spring-mass system, interacting via collisions only with a smooth particle system.