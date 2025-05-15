Fluid sim implemented using Stam's Real time fluid sim algorithm
http://graphics.cs.cmu.edu/nsp/course/15-464/Fall09/papers/StamFluidforGames.pdf

Optimized to use AVX2 intructions for maximum single thread performance.

Future Ideas
- [ ] Extend to 3D
- [ ] MultiThread
- [ ] Offload to GPU using compute shaders
- [ ] Fix Rendering to be a single texture being written to by the state, instead of NxN squares
