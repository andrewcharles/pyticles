A Force alters the rates of change of a ParticleSystem. It is the mapping from the particle attributes to the x and xdot state vectors that tells you which variables are fair game for being modified by a Force, and which aren't.

A Force is the right hand side of the differential equation xdot = f(x).

If a force depends on a rate of change then we have a potential ordering problem where the application of one force before another alters the outcome of the other force calculation. This is inconsistent, as we assume all the forces act simultaneously.

Most applications will only require the update of dv/dt and de/dt, and most forces will not have terms that depend on these values, so each force will just add to a force accumulator. It is probably not fruitful to try to design a more general system, capable of handling forces that depend on these rates of change (e.g. dv/dt = f(de/dt)).

A collision operator is a special class of force, because it is expected to alter the particles in place. When a collision is detected, velocities are changed immediately, and positions are adjusted for interpenetration. Because of this, a collision will not be a kind of force, but rather a different class of controller.