# Tiny-MCAB
##### <i>A smol version of tic-tac-toe MCAB with little-to-no optimization.</i>
Monte Carlo Alpha Beta is a search algorithm for making decisions in two-player
zero-sum games. 
# Algorithm Description
MCAB is an iterative algorithm that executes four subroutines until time runs out.

```
1. Selection
2. Expansion
3. Simulation
4. Back-propagation
```

## Selection
Starting at the root node, MCAB navigates the search tree in memory until reaching a node on the fringe.
At each internal node, a child is chosen according to some "tree policy." The tree policy should take into 
account the probability of winning stored in each child. A common policy is UCT, which treats each node 
like a multi-armed-bandit, exploring less and exploiting more as nodes accrue visits.

## Expansion
Once a leaf node has been selected, MCAB expands it into its children. A new node is added to the tree for
each legal action according to the leaf's state.

## Simulation
MCAB runs a shallow alpha-beta search below each new node to collect winning probabilities.
- 0   means that a player loses.
- 0.5 means that a player draws.
- 1   means that a player wins.
The probability may be any value in-between

## Back-propagation
Once probabilities have been determined, MCAB propagates them up the tree, following the unique path from 
the leaf to the root. The winning probability is added to the value of each node owned by the relevant 
alliance. The simulation count is added to the visit count of each node on the path, regardless of alliance.



