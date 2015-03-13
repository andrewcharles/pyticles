## Particle Refinement and Amalgamation ##
Suggest not allowing refinement every step - otherwise we can't study the effects of amalgamation shock. Should be a user control to initiate a refinement event. For long runs this could just be converted to a periodical refinement, or even monte-carlo refinement.

## Integration with existing Fortran modules ##
For a while a pure python rewrite was tempting. However this would be of little use
for anything but a toy implementation, unless I managed some neat numpy tricks which would destroy the readibility and ease of use gained from using OO structures in the first place.

As of now I have wrapped, and apparently sucessfully integrated the fortran 3d collisions module. It's a good feeling. Not only will this allow me to swap versions of components, and do cross testing for those components that are duplication, it will keep the existing Fortran code fresh and maintained.