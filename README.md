# Memory-integration
This is the simulation code used in the paper: "[Retrospective memory integration is accompanied with dynamic organization of cell assemblies]"

## Set up
You need Matlab to run this code. We tested it on Matlab R2018b.

## Usage
In the following we show how to reproduce the results reported in the paper. Raw simulation data will be stored in the D:\Memory_integration\data\.

1. Run gene_orthogonal_input_coding_neuron.m to generate orthogonal coding pattern for input layer.
2. Run gene_stim_patterns.m to generate stimuli file.
3. Run main.m to get the results of firing pattern, spike train and weight matrix.
4. Run weight_and_overlap_analysis.m to get the weight and coding neural overlap between pattern pairs of different degree across single item presentation.
5. Run Plot_weight_change.m to plot the weight change of pattern pairs of different degree across single item show time.
6. Run Plot_overlap_change.m to plot the neural coding overlap change of pattern pairs of different degree across single item show time.
7. Run Plot_firing_pattern_and_spike_train.m to get the firing pattern, spike train, inter-spike interval (ISI) and coeffecient of variation (CV) under different phases.

## References

Ye Wang et al., Retrospective memory integration is accompanied with dynamic organization of cell assemblies
