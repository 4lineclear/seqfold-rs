# callchain

This file is for listing how the recursion within this library works.

## Functions:

- fold
  - cache
  - traceback
- dg
  - fold
- dg_cache
  - cache
- dot_bracket
- cache
  - w
- w
  - w
  - v
  - min_value
  - multi_branch
- v
  - v
  - calc_pair
  - min_value
  - calc_stack
  - hairpin
  - bulge
  - internal_loop
  - multi_branch
- calc_pair
- min_value
- calc_d_g
- calc_j_s
- calc_stack
  - calc_pair
  - calc_d_g
- hairpin
  - calc_pair
  - calc_d_g
  - calc_j_s
- bulge
  - calc_pair
  - calc_d_g
  - calc_j_s
  - calc_stack
- internal_loop
  - calc_pair
  - calc_d_g
  - calc_j_s
  - calc_stack
- multi_branch
  - w
  - calc_stack
  - add_branch
- add_branch
  - w
  - add_branch
- traceback
  - traceback
  - trackback_energy
- trackback_energy

## Loops

Self-Loops

- w -> w
- v -> v
- add_branch -> add_branch
- traceback -> traceback

Cycles

- w -> multibranch -> w
- w -> v -> multibranch -> w
- w -> v -> multibranch -> add_branch -> w
