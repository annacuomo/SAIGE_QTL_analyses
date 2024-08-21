# Other tools comparison

## Direct comparisons

* [Nathan et al]() use a Poisson mixed model as implemented in R (add details).
* [CellRegMap]() has an `association test` mode that uses a linear mixed model as implemented in the [limix python package]().

Because these test for standard eQTL effects, we can compare the performance of SAIGE-QTL to these directly, both in terms of computation and memory cost, and to compare the actual results.

## Dynamic eQTL methods

On the other hand, [GASPACHO]() and CellRegMap's main test (the `interaction test`) are testing for something different, namely dynamic / (continuous) context-specific eQTL effects.
As the authors of both manuscripts state in their discussion, these tools do not scale well to large datasets, and should mostly be used for interpretation of existing eQTLs rather than for discovery.
More to the point, they test for different things, not the presence or absence of an eQTL but rather whether an eQTL's strength is modulated by a cell context / state.

We ran GASPACHO and CellRegMap's interaction test for OneK1K using PEER factors as cell states and show their computing time but a direct comparison with SAIGE-QTL would not be meanigful.
